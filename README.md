# Differential assembly diversifies GABA-A receptor structures and signaling
Scripts for the manuscript "Differential assembly diversifies GABA-A receptor structures and signaling"

Data is available at ftp://ftp.mrc-lmb.cam.ac.uk/pub/asente/differential-assembly/

Columns 1-3 in the file simulation_aa9_all.txt represent relative abundances for the three subunits (arbitrarily denoted as alpha/beta/x, where x can be gamma, delta or any other subunit).
Columns 4-12 represent relative interface likelihoods in the order: [aa, ab, ad, ba, bb, bd, da, db, dd] as given by the following matrix:

 Interface matrix:
    a b  d\n
 a aa ab ad\n
 b ba bb bd\n
 d da db dd\n
 
|   |  a  | b   | x  |
| --| --- | --- | -- |
| a | aa  | ab  | ad |
| b | ba  | bb  | bd |
| x | da  | db  | dd |

 First subunit in the pairwise interface is principal, second is complementary [ab = a+/b-]
 All receptor arrangements are read counterclockwise.

Columns 13-63 represent the abundances of the 51 receptor subtypes that can assemble from a pool three diffeent subunits. A lookup dataframe is generated within the analysis script.
