Function galecomplexity is given three integer parameters: 
 - e > 1 : dimension of the Gale diagram.
 - n > 0 : number of vectors on the Gale diagram.
 - m >= 0 : maximum value of the coordinates.

First, all possible e-dimensional Gale diagrams with n points and coordinates whose
absolute values are <= m are computed. Secondly, it is  checks that configurations
are full dimensional. Thirdly, configurations with interior points are discarded. 
It is done by checking that no hyperplane spanned by e-1 of the vectors strictly 
separates exactly one vector from the others. Finally, before adding a new 
configuration, it is checked that it has not the same face lattice as any other
already stored.

**Instructions

The code needs polymake instaled. To compile it, place the code in a folder F. Create
there a folder called "code" for example. In the terminal and from F execute polymake.
Then, *inside* the polymake shell (i.e., not from the command line), say
>found_extension("code");
>extend_application("/your/full/path/to/code");
Then in ../F/code/apps/polytope/src copy the code (the gale.cc file). This is the file
that should be modified in case you want to. A Makefile is then created by polymake.
To compile, in the Terminal and from ../F/code say
>make

Now inside the polymake shell you can calculate the configurations with
>gale_complexity(e,n,m);

If you only want the configurations belonging those with parameters (e,n,m) but not
with (e,n,m-1), you can calculate them with
>difference(e,n,m)

