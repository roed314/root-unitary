
 LABEL????
*isog_label                   | text      | label for isogeny class
*frac_ideal                   | numeric[] | basis for fractional ideal, expressed in terms of V^{g-1},...,V,1,F,F^2,...,F^g
*rep_type                     | smallint  | 0=ordinary or Centeleghe-Stix,...
*is_reduced                   | boolean   | Whether the fractional ideal is reduced (HNF, minimal norm, lexicographic within same norm)
 cm_type                      | boolean[] | Whether the +imaginary embedding is a p-adic non-unit, for embeddings sorted by real part
*cm_elt                       | numeric[] | An element of Q[F] that is positive imaginary under each embedding in the CM type
 can_be_principally_polarized | boolean   | Whether this abelian variety has a principal polarization
 rational_invariants          | numeric[] | Invariant factors of A(F_q)
*is_product                   | boolean   | Whether this isomorphism class is a product of smaller dimensional abelian varieties
 product_factorization        | jsonb     | List of pairs (label, e) expressing this as a product of smaller dimensional abelian varities (NULL if not)
 endo_ring                    | jsonb     | Some kind of description....
 related_objects              | text[]    | List of URLs



LABEL????
isom_label          | text       |
degree              | smallint   | degree of the polarization
kernel              | smallint[] | invariant factors for the kernel of the isogeny (cokernel of the map of lattices)
is_decomposible     | boolean    | Whether this polarized abelian variety is a product
decomposition       | jsonb      | List of pairs (label, e) expressing this polarized abelian variety as a product (NULL if not)
aut_group           | text       | GAP id
geom_aut_group      | text       | GAP id
is_serre_obstructed | smallint   | -1 if not a Jacobian, 0 if a hyperelliptic Jacobian, 1 if a nonhyperelliptic Jacobian
invariants          | jsonb      | For small genus, a list of geometric invariants (e.g. Igusa).  Only possible in the principal case







Things to add for the isogeny class

order_is_bass    | boolean | whether all the over-orders for the order Z[F,V] are Gorenstein
order_is_maximal | boolean | whether the order Z[F,V] is maximal
size             | integer | number of isomorphism classes within the isogeny class


For each isogeny class, write lines to two files
isomorphism_classes.txt (one line per ideal)
isog_label:frac_ideal:rep_type::cm_elt:is_product
e.g.
1.251.v:{{1,0},{0,1}}:0:f:{21,2}:f

isogeny_classes.txt (one line per class)
isog_label:order_is_bass:order_is_maximal:size
1.251.v:t:t:9

* Whether or not two isogeny classes come together after base extension, and what the degree is
* isogeny_graphs
* ideal class monoid (as an integer matrix)

