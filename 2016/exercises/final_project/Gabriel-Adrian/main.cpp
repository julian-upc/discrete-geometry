// Enumeration of simplicial polytopes and neighbourly
// simplicial polytopes in R⁴ with 8 vertices.
//
// Authors: Gabriel Riera Roca and Adrián Ranea Robles.
//
// Description:
//   This program lists 37 simplicial polytopes and 3
//   neighbourly simplicial polytopes in R⁴ with 8 vertices
//   by generating random affine Gale diagrams.
//
//   For example, with the seed 1234567890, approximately
//   32000 random diagrams are needed to list 3 neighbourly
//   simplicial polytopes and 37 simplicial polytopes.
//   This takes a couple of minutes in a personal computer.
//
//   Our list of polytopes agrees with the list obtained by
//   Grünbaum in his paper “An enumeration of simplicial
//   4-polytopes with 8 vertices”.
//
// Dependencies:
//   GNU Multiple Precision Arithmetic Library
//   Parma Polyhedra Library
//
#include <cstdint>
#include <array>
#include <vector>

#include <gmpxx.h>
#include <ppl.hh>
namespace ppl = Parma_Polyhedra_Library;

constexpr auto MAX_ITER = 1000000;

// Basis for R²
const auto X = ppl::Variable(0);
const auto Y = ppl::Variable(1);

// A point in R²
struct Point
{
    mpq_class x, y;

    Point(mpq_class x, mpq_class y) : x(x), y(y) {}

    // Return a generator to be used with PPL (convenience)
    ppl::Generator generator() const
    {
        auto nx = x.get_num();
        auto dx = x.get_den();
        auto ny = y.get_num();
        auto dy = y.get_den();

        mpz_class z = dx*dy;
        return ppl::Generator::point(nx*dy*X + ny*dx*Y, z);
    }
};

// A point in a Gale diagram
struct GalePt
{
    enum Sign { Positive, Negative };

    Point pt;
    Sign sign;
};

// Representation of a Gale diagram
struct GaleDiagram
{
    std::vector<GalePt> points;
};

// Stores the facets of a polytope.
struct CombinatorialType
{
    // Compact representation of a subset of a given set of size tsize
    // tsize is never larger than 8 in our case
    struct PointSet
    {
        using Bits = std::uint32_t;
        Bits bits;
        size_t tsize; // total size

        // Check for inclusion
        bool contained_in(CombinatorialType::PointSet other) const
        {
            return (bits & other.bits) == bits;
        }

        // Test for membership of an element
        bool contains(size_t i) const
        {
            return bits & (1 << i);
        }

        // Lexicographical ordering for sorting, etc
        struct TotalOrder
        {
            bool operator()(PointSet a, PointSet b) const { return a.bits < b.bits; }
        };
    };

    using VertId = size_t;
    using Facet = PointSet;

    size_t n_vert;
    std::vector<Facet> facets; // sorted by TotalOrder
};

bool operator==(CombinatorialType::PointSet a, CombinatorialType::PointSet b)
{
    return a.tsize == b.tsize && a.bits == b.bits;
}
bool operator!=(CombinatorialType::PointSet a, CombinatorialType::PointSet b)
{
    return !(a == b);
}


// Obtain the line through two points
ppl::Linear_Expression mk_line(Point p, Point q)
{
    mpq_class x3 = p.y - q.y;
    mpq_class y3 = q.x - p.x;
    mpq_class z3 = p.x*q.y - p.y*q.x;

    mpz_class d1 = x3.get_den();
    mpz_class d2 = y3.get_den();
    mpz_class d3 = z3.get_den();

    return d2*d3*x3.get_num()*X + d1*d3*y3.get_num()*Y + d1*d2*z3.get_num();
}

// Return whether the diagram is a affine Gale diagram of a polytope.
template <bool Neighborly>
bool is_gale_diagram_of_polytope(const GaleDiagram& diagram)
{
    const auto& points = diagram.points;

    for (size_t i = 0; i < points.size(); i++) {
        auto p = points[i].pt;
        for (size_t j = i+1; j < points.size(); j++) {
            auto q = points[j].pt;

            auto l = mk_line(p, q);

            // only used if Neighborly
            ppl::Generator_System positive;
            ppl::Generator_System negative;

            size_t n_left[2] = {0, 0};
            size_t n_right[2] = {0, 0};

            for (size_t k = 0; k < points.size(); k++) {
                if (k == i || k == j)
                    continue;

                auto r = points[k].pt;
                auto sign = points[k].sign;

                mpq_class scalar_product = r.x*l.coefficient(X) + r.y*l.coefficient(Y) + l.inhomogeneous_term();

                if (scalar_product < 0)
                    n_left[sign]++;
                else if (scalar_product > 0)
                    n_right[sign]++;

                if (Neighborly) {
                    if (sign == GalePt::Positive)
                        positive.insert(r.generator());
                    else
                        negative.insert(r.generator());
                }
            }

            if (n_left[GalePt::Positive] + n_right[GalePt::Negative] < 2)
                return false;

            if (n_left[GalePt::Negative] + n_right[GalePt::Positive] < 2)
                return false;

            if (Neighborly) {
                if (positive.empty() || negative.empty())
                    return false;

                auto tmp = ppl::C_Polyhedron(positive);
                tmp.intersection_assign(ppl::C_Polyhedron(negative));
                if (tmp.is_empty())
                    return false;
            }
        }
    }

    return true;
}

// Checks if a set of points is contained in a facet of the polytope
// represented by the diagram.
bool is_valid_face_candidate(const GaleDiagram& diagram, CombinatorialType::PointSet points)
{
    ppl::Generator_System positive;
    ppl::Generator_System negative;

    for (size_t i = 0; i < diagram.points.size(); i++) {
        if (!points.contains(i)) {
            if (diagram.points[i].sign == GalePt::Positive)
                positive.insert(diagram.points[i].pt.generator());
            else
                negative.insert(diagram.points[i].pt.generator());
        }
    }

    if (positive.empty() || negative.empty())
        return false;

    auto tmp = ppl::C_Polyhedron(positive);
    tmp.intersection_assign(ppl::C_Polyhedron(negative));
    return !tmp.is_empty();
}

// Obtains the combinatorial type of the polytope represented by the given
// diagram.
CombinatorialType get_combinatorial_type(const GaleDiagram& diagram)
{
    CombinatorialType result;
    result.n_vert = diagram.points.size();

    std::vector<CombinatorialType::PointSet> candidates;

    assert(diagram.points.size() < 8*sizeof(size_t));
    const size_t powerset = 1 << diagram.points.size();
    candidates.reserve(powerset);

    // powerset-1 because we don't want the whole polygon
    for (CombinatorialType::PointSet::Bits bits = 0; bits < powerset-1; bits++) {
        auto candidate = CombinatorialType::PointSet {
            .bits = bits,
            .tsize = diagram.points.size()
        };

        if (is_valid_face_candidate(diagram, candidate))
            candidates.push_back(candidate);
    }

    auto implies = [](bool p, bool q) { return !p || q; };
    result.facets.reserve(candidates.size()/2);

    for (auto f : candidates) {
        auto maximal =
            std::all_of(candidates.begin(), candidates.end(), [=](auto c) {
                return implies(f.contained_in(c), f == c);
            });

        if (maximal)
            result.facets.push_back(f);
    }

    sort(result.facets.begin(), result.facets.end(),
        CombinatorialType::PointSet::TotalOrder());
    return result;
}

// Checks whether two combinatorial types are equal by testing all relabeling
// of vertices in the facets of ct2.
bool are_isomorphic(const CombinatorialType& ct1, const CombinatorialType& ct2)
{
    if (ct1.n_vert != ct2.n_vert || ct1.facets.size() != ct2.facets.size())
        return false;

    auto n = ct1.n_vert;

    std::vector<CombinatorialType::VertId> perm;
    perm.resize(n);

    for (size_t i = 0; i < n; i++)
        perm[i] = i;

    CombinatorialType ct2_perm = ct2;

    do {
        for (size_t i = 0; i < ct2.facets.size(); i++) {
            const auto& ct2_facet = ct2.facets[i];

            auto& ct2_facet_perm = ct2_perm.facets[i];
            ct2_facet_perm.tsize = ct2_facet.tsize;
            ct2_facet_perm.bits = 0;

            for (size_t j = 0; j < n; j++)
                ct2_facet_perm.bits |= (size_t)ct2_facet.contains(j) << perm[j];
        }

        std::sort(ct2_perm.facets.begin(), ct2_perm.facets.end(),
                CombinatorialType::PointSet::TotalOrder());

        if (ct1.facets == ct2_perm.facets)
            return true;
    } while (std::next_permutation(perm.begin(), perm.end()));

    return false;
}

// Return whether it is the combinatorial type of a simplicial polytope
bool is_simplicial_polytope(const CombinatorialType& ct, size_t dim)
{
    return std::all_of(ct.facets.begin(), ct.facets.end(), [dim](auto f) {
        auto bits = f.bits;
        size_t count = 0;

        do {
            count += bits & 1;
        } while (bits >>= 1);

        return count == dim;
    });
}


// Return the generator of random integers.
gmp_randclass& get_random_gen()
{
    static gmp_randclass* rnd = [] {
        static gmp_randclass rnd (gmp_randinit_default);
        rnd.seed(1234567890);
        return &rnd;
    }();

    return *rnd;
}

// Return a random affine GaleDiagram whose points are
// quotient of random integers of size bits.
template <size_t N>
GaleDiagram get_random_gale_diagram(mpz_class bits=32)
{
    auto& rnd = get_random_gen();
    mpz_class rndsigns = rnd.get_z_bits(N);

    GaleDiagram diagram;

    for (size_t i = 0; i < N; i++) {
        mpq_class x = {rnd.get_z_bits(bits), rnd.get_z_bits(bits)};
        mpq_class y = {rnd.get_z_bits(bits), rnd.get_z_bits(bits)};

        // Test bit bit_index in op and return 0 or 1 accordingly.
        auto sign =
            mpz_tstbit(rndsigns.get_mpz_t(), i)
                ? GalePt::Positive
                : GalePt::Negative;

        diagram.points.push_back(GalePt {
            .pt = Point(x, y),
            .sign = sign
        });
    }

    return diagram;
}

// Iterates over randomly generated affine Gale diagrams of N points.
// Produces MAX_ITER diagrams.
template <size_t N, typename Yield>
void foreach_gale_diagram(Yield&& yield)
{
    for (size_t i = 0; i < MAX_ITER; i++) {
        const auto diagram = get_random_gale_diagram<N>();
        yield(diagram);

        if (i > 0 && i % (MAX_ITER/10) == 0)
            std::cout << "( " << 100*i/MAX_ITER << "% )\n";
    }
}

// Iterates over (unique) combinatorial types of (N-4)-polytopes of N vertices.
// If Neighborly is true, only neighborly polytopes are considered.
template <size_t N, bool Neighborly, typename Yield>
void foreach_combinatorial_type(Yield&& yield)
{
    std::vector<CombinatorialType> types;

    foreach_gale_diagram<N>([&](auto&& diagram) {
        if (!is_gale_diagram_of_polytope<Neighborly>(diagram))
            return;

        const auto new_type = get_combinatorial_type(diagram);

        if (!is_simplicial_polytope(new_type, N-4))
            return;

        auto repeated =
            std::any_of(types.begin(), types.end(), [&](auto&& type) {
                return are_isomorphic(type, new_type);
            });

        if (!repeated) {
            types.push_back(new_type);
            yield(diagram, new_type);
        }
    });
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////   TESTS BELOW   /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


// Print a combinatorial type in human-readable format (list of facets)
void print_combinatorial_type(const CombinatorialType& ct)
{
    std::cout << "# facets: " << ct.facets.size() << "\n";

    for (auto facet : ct.facets) {
        for (size_t i = 0; i < facet.tsize; ++i)
            if (facet.contains(i))
                std::cout << i;
        std::cout << " ";
    }

    std::cout << std::dec << std::endl;
}

// The Gale diagram of C4(8) for testing purposes
GaleDiagram get_c48_gale_diagram()
{
    return GaleDiagram {
        .points = {
            { .pt = Point( 2,  1), .sign = GalePt::Negative },
            { .pt = Point( 1,  2), .sign = GalePt::Positive },
            { .pt = Point(-1,  2), .sign = GalePt::Negative },
            { .pt = Point(-2,  1), .sign = GalePt::Positive },
            { .pt = Point(-2, -1), .sign = GalePt::Negative },
            { .pt = Point(-1, -2), .sign = GalePt::Positive },
            { .pt = Point( 1, -2), .sign = GalePt::Negative },
            { .pt = Point( 2, -1), .sign = GalePt::Positive }
        }
    };
}

// The counter example in the Ziegler book, exercise 6.15 for testing purposes
GaleDiagram get_non_c48_gale_diagram()
{
    return GaleDiagram {
        .points = {
            { .pt = Point( 2,  3), .sign = GalePt::Negative },
            { .pt = Point( 1,  2), .sign = GalePt::Positive },
            { .pt = Point(-1,  2), .sign = GalePt::Negative },
            { .pt = Point(-2,  3), .sign = GalePt::Positive },
            { .pt = Point(-2, -1), .sign = GalePt::Negative },
            { .pt = Point(-1, -2), .sign = GalePt::Positive },
            { .pt = Point( 1, -2), .sign = GalePt::Negative },
            { .pt = Point( 2, -1), .sign = GalePt::Positive }
        }
    };
}

// The combinatorial type of C4(8) for testing purposes
CombinatorialType get_c48_combinatorial_type()
{
    auto ct = CombinatorialType {
        .n_vert = 8,
        .facets = {
            { .bits = 0b00001111, .tsize = 8 }, // 0123
            { .bits = 0b00011110, .tsize = 8 }, // 1234
            { .bits = 0b00011011, .tsize = 8 }, // 0135
            { .bits = 0b00110011, .tsize = 8 }, // 0145
            { .bits = 0b00110110, .tsize = 8 }, // 1245
            { .bits = 0b00111100, .tsize = 8 }, // 2345
            { .bits = 0b01111000, .tsize = 8 }, // 3456
            { .bits = 0b01101100, .tsize = 8 }, // 2356
            { .bits = 0b01100110, .tsize = 8 }, // 1256
            { .bits = 0b01100011, .tsize = 8 }, // 0156
            { .bits = 0b11000011, .tsize = 8 }, // 0167
            { .bits = 0b11100001, .tsize = 8 }, // 0567
            { .bits = 0b11000110, .tsize = 8 }, // 1267
            { .bits = 0b11001100, .tsize = 8 }, // 2367
            { .bits = 0b10001101, .tsize = 8 }, // 0237
            { .bits = 0b10000111, .tsize = 8 }, // 0127
            { .bits = 0b10011001, .tsize = 8 }, // 0347
            { .bits = 0b10110001, .tsize = 8 }, // 0457
            { .bits = 0b11011000, .tsize = 8 }, // 3467
            { .bits = 0b11110000, .tsize = 8 }, // 4567
        }
    };

    std::sort(ct.facets.begin(), ct.facets.end(),
            CombinatorialType::PointSet::TotalOrder());

    return ct;
}

int main()
{
    // c48 is isomorphic to c48
    assert(
        are_isomorphic(
            get_combinatorial_type(get_c48_gale_diagram()),
            get_c48_combinatorial_type()));

    // c48 is not isomorphic to non-c48
    assert(
        !are_isomorphic(
            get_combinatorial_type(get_c48_gale_diagram()),
            get_combinatorial_type(get_non_c48_gale_diagram())));

    size_t cnt = 0;

    // change true to false to list all simplicial polytopes
    foreach_combinatorial_type<8, true>([&](auto&&, auto&& comb_type) {
        print_combinatorial_type(comb_type);
        cnt++;
    });

    std::cout << "\ntotal count: " << cnt << std::endl;
}
