./filter_vdjdb.R "HomoSapiens" "TRB" "MHCI" 1 10
mvn clean install
java -Xmx10G -jar target/cdr3align-0.0.1.jar 4 2 vdjdb.txt mutations.txt
cat mutations.txt | wc -l
rm mutations.txt.gz
gzip mutations.txt
