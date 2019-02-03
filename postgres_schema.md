Table name: `av_fqisog`

This table represents unpolarized abelian varieties, up to isogeny.

Things to add for the isogeny class

Column                       | Type       | Notes
-----------------------------|------------|------
label                        | text       |
g                            | smallint   |
q                            | integer    |
poly_str                     | text       | Space separated string of coefficients, for searching
p_rank                       | smallint   |
dim1_factors                 | smallint   | Number of dimension 1 factors
dim2_factors                 | smallint   | Number of dimension 2 factors
dim3_factors                 | smallint   | Number of dimension 3 factors
dim4_factors                 | smallint   | Number of dimension 4 factors
dim5_factors                 | smallint   | Number of dimension 5 factors
dim1_distinct                | smallint   | Number of distinct dimension 1 factors
dim2_distinct                | smallint   | Number of distinct dimension 2 factors
dim3_distinct                | smallint   | Number of distinct dimension 3 factors
dim4_distinct                | smallint   | Number of distinct dimension 4 factors
dim5_distinct                | smallint   | Number of distinct dimension 5 factors
poly                         | integer[]  | Coefficients of the Weil polynomial.  The first will always be `1` and the last `q^g`
angles                       | float8[]   | Angles corresponding to roots in the closure of the upper half plane, divided by `pi`.  All will be in the interval `[0, 1]`, and there will be `g` of them unless `0` or `1` is included.
angle_rank                   | smallint   | The dimension of the Q-span of the angles (see knowl for complete definition)
slopes                       | text[]     | The sorted list of slopes, as string representations of rational numbers.  Duplicated slopes will have "A", "B", etc appended.
abvar_counts                 | numeric[]  | The list of counts `#A(F_{q^i})` for `i=1..10`, for `A` in this isogeny class
abvar_counts_str             | text       | A space separated string of abelian variety counts, for searching
curve_counts                 | numeric[]  | The list of curve counts `#C(F_{q^i})` for `i=1..10` for any curve `C` of genus `g` with `J(C)` in this isogeny class
curve_counts_str             | text       | A space separated string of curve counts, for searching
point_count                  | integer    | The count `#C(F_q)`, duplicated for searching purposes
has_jacobian                 | smallint   | 1 if it is known that this isogeny class contains a Jacobian; -1 if it is known that it does not; 0 otherwise
has_principal_polarization   | smallint   | 1 if it is known that this isogeny class contains a principally polarizable abelian variety; -1 if it is known that it does not; 0 otherwise
is_simple                    | boolean    |
simple_factors               | text[]     | A list of labels of simple factors.  Duplicated factors will have "A", "B", etc appended.
simple_distinct              | text[]     | A list of distinct labels of simple factors.
simple_multiplicities        | smallint[] | For each distinct simple factor, the multiplicity in the decomposition.
number_field_degrees         | smallint[] | For each distinct simple factor, the degree of the the corresponding number field.
divalg_dimensions            | smallint[] | For each distinct simple factor, the dimension of the division algebra over the number field.
number_fields                | text[]     | The number fields associated to the irreducible factors of the Weil polynomial
galois_groups                | text[]     | The Galois groups of the number fields associated to the irreducible factors of the Weil polynomial, e.g. "4T3"
places                       | text[]     | A list of lists of lists of rational numbers stored as strings, giving the prime ideals above `p`.  The terms in the outer list correspond to distinct simple factors, the terms in the middle lists correspond to places in the corresponding number field, and each inner list gives coefficients for a two-element generator of that prime ideal (along with `p`) as coefficients of powers of `F`.
brauer_invariants            | text[]     | A list of lists of rational numbers stored as strings.  The terms in the outer list correspond to distinct simple factors, and the terms in each inner list correspond to the places in the corresponding number field.
geometric_extension_degree   | smallint   | The smallest degree extension of the base field over which the endomorphism algebra becomes the full endomorphism algebra
geometric_simple_factors     | text[]     | A list of labels of simple factors after base changing by the geometric extension degree.  Duplicated factors will have "A", "B", etc appended.  NULL if geometric_extension_degree is 1.
geometric_simple_distinct    | text[]     | A list of distinct labels of simple factors after geometric base change.  NULL if geometric_extension_degree is 1.
geometric_multiplicities     | smallint[] | For each distinct geometric simple factor, the multiplicity in the decomposition.  NULL if geometric_extension_degree is 1.
geometric_number_field_degrees| smallint[] | For each distinct geometric simple factor, the degree of the corresponding number field.  NULL if geometric_extension_degree is 1.
geometric_divalg_dimensions  | smallint[] | For each distinct geometric simple factor, the dimension of the division algebra over the number field.  NULL if geometric_extension_degree is 1.
geometric_number_fields      | text[]     | The number fields associated to the irreducible factors of the base-changed Weil polynomial.  NULL if geometric_extension_degree is 1.
geometric_galois_groups      | text[]     | The Galois groups of the geometric number fields.  NULL if geometric_extension_degree is 1.
geometric_places             | text[]     | Places, after base changing to the geometric field.
geometric_brauer_invariants  | text[]     | Brauer invariants, after bas changing to the geometric field.
primitive_models             | text[]     | A list of labels giving primitive models for this isogeny class (ie, this class arises from base change from the model).  If primitive, NULL.
is_primitive                 | boolean    |
twists                       | jsonb      | A list of triples `(label, geom_label, r)` where `label` is the label of a twist, `r` is an extension degree where the twists become isomorphic, and `geom_label` is the label of the common base change to that degree.
size                         | integer    | number of isomorphism classes within the isogeny class
zfv_is_bass                  | boolean    | whether all the over-orders for the order `Z[F,V]` are Gorenstein
zfv_is_maximal               | boolean    | whether the order `Z[F,V]` is maximal
zfv_index                    | numeric    | the index of the order `Z[F,V]` in the maximal order
zfv_index_factorization      | numeric[]  | A list of pairs (p, e) giving the factorization of the index
zfv_plus_index               | numeric    | the index of the order `Z[F+V]` in the maximal order of the real subfield
zfv_plus_index_factorization | numeric[]  | A list of pairs (p, e) giving the factorization of the index
zfv_plus_norm                | numeric    | The absolute value of the norm of F-V to Z
isogeny_graphs               | jsonb      | list of pairs `(p, G)`, where `p` is a degree (or maybe list of degrees) and `G` is a list of pairs `(u,v)` representing the directed edge from `u` to `v`.  Each of `u` and `v` is the `isom_letter` for the corresponding isogeny class
ideal_class_generators       | text[]     | A list of `isom_letters` for isomorphism classes that generate the ideal monoid
ideal_class_relations        | integer[]  | A matrix of positive integers giving relations between the ideal class generators



Table name: `av_fqisom`

This table represents unpolarized abelian varieties, up to isomorphism.

Column                       | Type      | Notes
-----------------------------|-----------|------
label                        | text      | `g.q.weil.enum`, where `g` is the dimension, `q` is the cardinality of the base field, `weil` is the encoding of the Weil polynomial and `enum` is `isom_letter`
isom_num                     | integer   | A 0-based enumeration of the isomorphism classes within an isogeny class, TBD
isom_letter                  | text      | Base 26 a-z encoding of isom_num
isog_label                   | text      | label for isogeny class
frac_ideal_numerators        | numeric[] | numerators for a basis for the fractional ideal, expressed in terms of V^{g-1},...,V,1,F,F^2,...,F^g
frac_ideal_denominators      | numeric[] | denominators for a basis for the fractional ideal, expressed in terms of V^{g-1},...,V,1,F,F^2,...,F^g.  NULL if all 1.
over_order                   | text      | isom_letter for over order associated to this ideal MAYBE???
is_over_order                | boolean   | MAYBE???
weak_equivalence_class       | text      | label for the weak equivalence class???
rep_type                     | smallint  | 0=ordinary or Centeleghe-Stix,...
is_reduced                   | boolean   | Whether the fractional ideal is reduced (HNF, minimal norm, lexicographic within same norm)
cm_type                      | boolean[] | Whether the +imaginary embedding is a p-adic non-unit, for embeddings sorted by real part
cm_elt                       | numeric[] | An element of Q[F] that is positive imaginary under each embedding in the CM type
can_be_principally_polarized | boolean   | Whether this abelian variety has a principal polarization
rational_invariants          | numeric[] | Invariant factors of A(F_q)
is_product                   | boolean   | Whether this isomorphism class is a product of smaller dimensional abelian varieties
product_factorization        | jsonb     | List of pairs (label, e) expressing this as a product of smaller dimensional abelian varities (NULL if not)
endo_ring                    | jsonb     | Some kind of description....
related_objects              | text[]    | List of URLs


Table name: `av_fqpol`

This table represents polarized abelian varieties, up to isomorphism.

Column              | Type       | Notes
--------------------|------------|------
label               | text       | ?????
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

Column                       | Type      | Notes
-----------------------------|-----------|------
size                         | integer   | number of isomorphism classes within the isogeny class
zfv_is_bass                  | boolean   | whether all the over-orders for the order `Z[F,V]` are Gorenstein
zfv_is_maximal               | boolean   | whether the order `Z[F,V]` is maximal
zfv_index                    | numeric   | the index of the order `Z[F,V]` in the maximal order
zfv_index_factorization      | numeric[] | A list of pairs (p, e) giving the factorization of the index
zfv_plus_index               | numeric   | the index of the order `Z[F+V]` in the maximal order of the real subfield
zfv_plus_index_factorization | numeric[] | A list of pairs (p, e) giving the factorization of the index
zfv_plus_norm                | numeric   | The absolute value of the norm of F-V to Z
isogeny_graphs               | jsonb     | list of pairs `(p, G)`, where `p` is a degree (or maybe list of degrees) and `G` is a list of pairs `(u,v)` representing the directed edge from `u` to `v`.  Each of `u` and `v` is the `isom_letter` for the corresponding isogeny class
ideal_class_generators       | text[]    | A list of `isom_letters` for isomorphism classes that generate the ideal monoid
ideal_class_relations        | integer[] | A matrix of positive integers giving relations between the ideal class generators

 * Whether or not two isogeny classes come together after base extension, and what the degree is

* Change `brauer_invs` from a string to list of lists of strings
* Change jsonb types to arrays
* Change `nf` to `text[]`: a list of number fields in the non-simple case
* Update galois_t when not set
* Add missing isogeny classes (t^2-p)
* Rename/retype some columns:
 - poly (jsonb -> integer[])
 - angles (jsonb -> float8[])
 - ang_rank -> angle_rank (smallint)
 - slps -> slopes (jsonb -> text[])
 - A_cnts -> abvar_counts (jsonb -> numeric[])
 - A_cnts_str -> abvar_counts_str (text)
 - C_cnts -> curve_counts (jsonb -> numeric[])
 - C_cnts_str -> curve_counts_str (text)
 - pt_cnt -> point_count (integer)
 - is_jac -> has_jacobian (smallint -> boolean)
 - is_pp -> has_principal_polarization (smallint -> boolean)
 - decomp -> decomposition (jsonb)
 - is_simp -> is_simple (boolean)
 - simple_factors (jsonb -> text[])
 - simple_distinct (jsonb -> text[])
 - brauer_invs -> brauer_invariants (text -> text[])
 - places (jsonb -> text[])
 - prim_models -> primitive_models (jsonb -> text[])
 - is_prim -> is_primitive (boolean)
 - nf -> number_fields (text -> text[])
 - galois_t -> galois_groups (smallint -> text[])
 - galois_n -> XXX

Unchanged:
 - label (text)
 - g (smallint)
 - q (integer)
 - poly_str (text)
 - p_rank (smallint)
 - dim1_factors (smallint)
 - dim2_factors (smallint)
 - dim3_factors (smallint)
 - dim4_factors (smallint)
 - dim5_factors (smallint)
 - dim1_distinct (smallint)
 - dim2_distinct (smallint)
 - dim3_distinct (smallint)
 - dim4_distinct (smallint)
 - dim5_distinct (smallint)

June Ju and Everett Howe: telling whether an abelian variety is absolutely simple
Look at primes dividing discriminant of Weil field, products at most 4g^2.
Find pairwise r, hash on multiple r
Supersingular if and only if the ultimate field is Q

For Stefano
-----------

We'll be doing isogeny classes that Stefano has already computed, plus:

 * Any g, q
 * ordinary or C-S (q=p, no real roots), squarefree
 * Z[F,V] = maximal order

For each isogeny class, write lines to two files

 * isomorphism_classes.txt (one line per ideal)
```
isog_label:frac_ideal:rep_type:is_reduced:cm_elt:is_product
e.g.
1.251.v:{{1,0},{0,1}}:0:f:{21,2}:f
```
 * isogeny_classes.txt (one line per class)
```
isog_label:order_index:order_is_bass:order_is_maximal:size
e.g.
1.251.v:1:t:t:9
```

Use `\N` for null.

["1.251.v", [251,21,1]]