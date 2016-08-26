# uclust --input /tmp/UclustExactMatchFilterfxW19_.fasta --id 0.97 --tmpdir /tmp --optimal --w 8 --stepwords 8 --usersort --maxaccepts 1 --stable_sort --maxrejects 8 --uc data/otus_97_bysample/85_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	297	*	*	*	*	*	QiimeExactMatch.85_1483179	*
S	1	300	*	*	*	*	*	QiimeExactMatch.85_1483175	*
S	2	306	*	*	*	*	*	QiimeExactMatch.85_1483181	*
H	0	297	98.7	+	0	0	297M	QiimeExactMatch.85_1483176	QiimeExactMatch.85_1483179
H	0	293	99.3	+	0	0	91M2I31M2I171M	QiimeExactMatch.85_1483177	QiimeExactMatch.85_1483179
H	1	300	98.3	+	0	0	300M	QiimeExactMatch.85_1483173	QiimeExactMatch.85_1483175
H	1	300	98.3	+	0	0	300M	QiimeExactMatch.85_1483174	QiimeExactMatch.85_1483175
H	1	299	98.3	+	0	0	208MI91M	QiimeExactMatch.85_1483182	QiimeExactMatch.85_1483175
H	0	295	99.0	+	0	0	91M2I204M	QiimeExactMatch.85_1483178	QiimeExactMatch.85_1483179
H	0	293	99.7	+	0	0	91M2I31M2I171M	QiimeExactMatch.85_1483183	QiimeExactMatch.85_1483179
C	0	5	99.2	*	*	*	*	QiimeExactMatch.85_1483179	*
C	1	4	98.3	*	*	*	*	QiimeExactMatch.85_1483175	*
C	2	1	*	*	*	*	*	QiimeExactMatch.85_1483181	*
