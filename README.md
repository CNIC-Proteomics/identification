# Script "rels2pq_unique"

## Description

The script **rels2pq_unique.rel_input.py** selects one protein for each peptide in the cases there are redundances.

The decision is based on the sum of peptides from a contained protein. In this way, for a list of proteins in one peptide, we select the protein with the higher number of peptides in total.

In the case of draw, we apply the following order of decisions:

1. We apply a preferenced text from a regular expression.

2. If there is not decision yet and the database has been given, we provide the protein with the sorted sequence.

3. If there is not decision yet, we provide the protein by alphabetic order.


## Input

Text file with the following columns (tab-separated values):
[LPP], [FASTAProteinDescription], [Sequence], [Redundances]

## Output

File with the following columns (tab-separated values):
    [FASTAProteinDescription], [Sequence], [Tags]

The "Tags" columns descriptions:

- LPP from the current peptide. At the moment, it is equal than the number of scans that appears in the data set.
- LPQ from the current protein (Sum of )
- Number of the step decision:
    0, decision by higher value of LPQ
    1, decision by the regular expression parameter
    2, decision by the sorted sequence
    3, decision by the alphabetic order


## Battery tests

### Execution
```
cd test

./run_rel2pq_unique.rel_input.sh
```

### Parameters

- test x.1

It is not take into account any parameters. In the case of draw, we provide the sorted protein.

- test x.2

Taking into account the species. Then, in the case of draw, we give preferences the proteins for the given species.

- test x.3

Taking into account a regular expression. For example, with ">sp" give a preferences to the SwissProt proteins.

- test x.4

Taking into account the lemgth of the sequences. Then, in the case of draw, we provide the sorted protein.

- test x.5

Taking into account the species, a regular expression, and the length of the sequences.

### Input datasets

The paths of following dataset are within the CNIC server (tierra.cnic.es).

#### 1. Mouse, DatosCrudos\Jose_Antonio_Enriquez\Ratones_Heteroplasmicos\NHF\LIVER

The ID-q and the relationship table, comes from the following path:
```
S:\U_Proteomica\UNIDAD\DatosCrudos\Jose_Antonio_Enriquez\Ratones_Heteroplasmicos\NHF\LIVER\Pre-SanXoT\ID-q.txt

echo "prepare working files"
cat header.tsv > p2q_rels_aux.tsv
awk -F "\t" '{print 1"\t"$9"\t"$10"\t"$20}' ID-q.txt | sed 1d >> p2q_rels_aux.tsv
```

The database comes from UniProt (2019_02) - organism:"Mus musculus (Mouse) [10090]" AND proteome:up000000589 (54,185 proteins)  -

<!-- Note: there were some mistakes in the replaces of modification. To fix this problem.
```
sed -e 's/K@/K*/g'  -e 's/K#/K*/g' -e 's/C\*/C@/g' -e 's/C#/C@/g' -e 's/M\*/M#/g' -e 's/M@/M#/g'    p2q_rels_aux.tsv > input_p2q_rels.tsv
```
-->









<!-- DEPRECATED -->


<!-- ### 1. Human, LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro

The ID-q and the relationship table, comes from the following path:
```
S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\CS_nitrations\SanXot_2\TMT1\ID-q-TMT1.txt
```
The database comes from:
```
S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\Databases\Homo_sapiens_release_swiss_prot_2018_07_Conc.fasta
``` -->

<!-- 2. LAB_JMR\Marfan\Mice\plasma (coming from COMET => 'Redundance' column does NOT EXIST)

The ID-q and the relationship table, comes from the following path:
```
S:\U_Proteomica\LABS\LAB_JMR\Marfan\Mice\plasma\Comet\DigPar\Sanxotp2P\Shifts-ID-Q_XsVs_vs_AveContr.txt
```
The database comes from:
```
S:\U_Proteomica\LABS\LAB_JMR\Alvaro\uniprot-Mus+musculus_release 2018_04_Crap_Tomato.fasta
```
 -->





<!-- # Check the list of peptides

cut -f 3 input_p2q_rels.tsv | sed 1d | sort -u | sed 's/"//g' > peptides_from_given_input.txt

cut -f 2 p2q_rels_unique.test1.tsv | sed 1d | sort -u > peptides_from_results.txt

comm -3 peptides_from_given_input.txt peptides_from_results.txt -->
