# uclust --input /tmp/UclustExactMatchFilterqV3QMG.fasta --id 0.97 --tmpdir /tmp --optimal --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc data/otus_97_bysample/84_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	306	*	*	*	*	*	QiimeExactMatch.84_44718	*
S	1	300	*	*	*	*	*	QiimeExactMatch.84_44719	*
H	1	300	98.7	+	0	0	300M	QiimeExactMatch.84_44721	QiimeExactMatch.84_44719
S	2	297	*	*	*	*	*	QiimeExactMatch.84_44722	*
H	2	295	99.3	+	0	0	91M2I204M	QiimeExactMatch.84_44720	QiimeExactMatch.84_44722
C	0	1	*	*	*	*	*	QiimeExactMatch.84_44718	*
C	1	2	98.7	*	*	*	*	QiimeExactMatch.84_44719	*
C	2	2	99.3	*	*	*	*	QiimeExactMatch.84_44722	*
