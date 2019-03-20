#!/usr/bin/bash

run_cmd() {
  echo $1
  time eval $1
}

INFILE="rel2pq_unique/input_p2q_rels.tsv"
OUTFILE="rel2pq_unique/p2q_rels_unique"
DBFILE="rel2pq_unique/uniprot-filtered-proteome_UP000000589+AND+organism__Mus+musculus+(Mouse%--.fasta"
#DEBUG=" -vv "

echo "** delete log files"
run_cmd "rm ${OUTFILE}.*.log"

# Test 1:
# It is not take into account any parameters
name="test1"
echo "** ${name}"
run_cmd "python ../rels2pq_unique.rel_input.py -i ${INFILE} -r ${OUTFILE}.${name}.tsv -l ${OUTFILE}.${name}.log ${DEBUG}"

# Test 2:
# Taking into account the species
name="test2"
echo "** ${name}"
run_cmd "python ../rels2pq_unique.rel_input.py -i ${INFILE} -r ${OUTFILE}.${name}.tsv -l ${OUTFILE}.${name}.log ${DEBUG} -s 'Mus musculus' "

# Test 3:
# Taking into account a regular expression
name="test3"
echo "** ${name}"
run_cmd "python ../rels2pq_unique.rel_input.py -i ${INFILE} -r ${OUTFILE}.${name}.tsv -l ${OUTFILE}.${name}.log ${DEBUG} -p '>sp'"

# Test 4:
# Taking into account the database file as parameter
name="test4"
echo "** ${name}"
run_cmd "python ../rels2pq_unique.rel_input.py -i ${INFILE} -r ${OUTFILE}.${name}.tsv -l ${OUTFILE}.${name}.log ${DEBUG} -d '${DBFILE}'"

# Test 5:
# Taking into account a regular expression
name="test5"
echo "** ${name}"
run_cmd "python ../rels2pq_unique.rel_input.py -i ${INFILE} -r ${OUTFILE}.${name}.tsv -l ${OUTFILE}.${name}.log ${DEBUG} -s 'Mus musculus' -p '>sp' -d '${DBFILE}'"


