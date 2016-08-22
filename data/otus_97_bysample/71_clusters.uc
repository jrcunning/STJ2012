# uclust --input /tmp/UclustExactMatchFilteraFRhX_.fasta --id 0.97 --tmpdir /tmp --optimal --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc data/otus_97_bysample/71_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	300	*	*	*	*	*	QiimeExactMatch.71_1271673	*
H	0	300	98.0	+	0	0	300M	QiimeExactMatch.71_1271674	QiimeExactMatch.71_1271673
H	0	300	98.3	+	0	0	300M	QiimeExactMatch.71_1271672	QiimeExactMatch.71_1271673
C	0	3	98.2	*	*	*	*	QiimeExactMatch.71_1271673	*
