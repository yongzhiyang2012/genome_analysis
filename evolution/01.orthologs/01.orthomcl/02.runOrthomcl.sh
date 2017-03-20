mysql -u root < mysql.setting.template
orthomclFilterFasta compliantFasta/ 10 20
makeblastdb -in goodProteins.fasta -dbtype prot
blastp -db goodProteins.fasta -query goodProteins.fasta -out all-all.blastp.out -evalue 1e-5 -outfmt 6 -num_threads 24
orthomclInstallSchema orthomcl.config.template
orthomclBlastParser all-all.blastp.out compliantFasta > similarSequences.txt
perl -p -i -e 's/0\t0/1\t-181/' similarSequences.txt
orthomclLoadBlast orthomcl.config.template similarSequences.txt
orthomclPairs orthomcl.config.template orthomcl_pairs.log cleanup=no
orthomclDumpPairsFiles orthomcl.config.template
mcl mclInput --abc -I 1.5 -o mclOutput
orthomclMclToGroups cluster 1 < mclOutput > groups.txt
