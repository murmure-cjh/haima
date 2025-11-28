import pandas as pd
import sys
import os
# usage: python deal_sma_workout.py /haplox/users/xiaoyue/16_sop/xwes/pipeline/SZ20240208064ORS-b/result/SMA/SZ20240208064ORS-b.smacacnv.csv SZ20240208064ORS-b.smacacnv.workout.csv SZ20240208064ORS-b.dedup.bam
# 1.input
sma_workout_inputf_wd=sys.argv[1] 
sma_workout_outputf_wd=sys.argv[2]
sample=sys.argv[3]
# 2.deal format
header = list(pd.read_csv(sma_workout_inputf_wd, sep=',', nrows=1, header=None).iloc[0])
data = pd.read_csv(sma_workout_inputf_wd, sep='|', skiprows=1,header=None,index_col=0)
#data.columns = header.iloc[0]
data.columns = header[1:len(header)+1]
data.index = [s.split('/')[-1] for s in list(data.index)]
# 3.judge
#===============================================================================================================
# Pi_p:	scaled proportion of SMN1 reads for positions p.
# cov_x_p:	raw coverage of gene x at position p.
# avg_cov_x:	average coverage for the whole gene x.
# std_control:	standard deviation for the average coverage of the 20 control genes.
# g.27134T>G:	consensus sequence at position 27134, as well as counts for "A", "C", "G" and "T".
# g.27706_27707delAT:	consensus sequence at positions 27706-27707, as well as counts for "A", "C", "G" and "T".
# scale_factor:	scale factor proportional to the total SMN1 and SMN2 copy-number.
# SMA carriers with a single SMN1 copy are expected to have Pi_b values under 1/3. However, complex SMN reorganizations may lead to large differences between Pi_a, Pi_b and Pi_c. \
# These cases should be analyzed carefully.
# The scale_factor, which is proportional to the absolute number of SMN1 and SMN2 copies, and cov_x_p can be used to estimate the absolute SMN1:SMN2 copy-number as follows:
# genotype	scale_factor	cov_SMN1_p/cov_SMN2_p
# 1:3	1	1/3
# 1:2	0.75	1/2
# 1:1	0.5	1
# In order to detect the so-called silent carriers (i.e. individuals with two copies of SMN1 on one chromosome, but none on the other), the consensus sequence at the two locations should also be taken into account. Depending on the number of SMN2 copies, the expected scale_factor should be close to 0.75 (2:1) or 0.5 (2:0) and, in both cases, the scaled proportion of SMN1 reads Pi_p should be close to 1/2 in each position.
#===============================================================================================================
# 这段文字主要是关于判断SMN1基因的拷贝数的方法。其中涉及到的一些术语和概念如下：
# Pi_p: 位置p的SMN1读数的比例。
# cov_x_p: 基因x在位置p的原始覆盖度。
# avg_cov_x: 基因x的平均覆盖度。
# std_control: 20个对照基因平均覆盖度的标准差。
# g.27134T>G: 位置27134的一致序列，以及"A"、"C"、"G"和"T"的计数。
# g.27706_27707delAT: 位置27706-27707的一致序列，以及"A"、"C"、"G"和"T"的计数。
# scale_factor: 与SMN1和SMN2总拷贝数成比例的缩放因子。
# SMA携带者有一个SMN1拷贝的预期Pi_b值小于1/3。然而，复杂的SMN重组可能导致Pi_a、Pi_b和Pi_c之间存在较大差异。这些情况应该仔细分析。
# scale_factor与SMN1和SMN2的绝对拷贝数成比例，cov_x_p可以用来估计绝对的SMN1:SMN2拷贝数，方法如下：
# 基因型：scale_factor：cov_SMN1_p/cov_SMN2_p
# 1:3：1：1/3
# 1:2：0.75：1/2
# 1:1：0.5：1
# 为了检测所谓的“静默携带者”（即在一条染色体上有两个SMN1拷贝，而在另一条染色体上没有），还应考虑两个位置的一致序列。根据SMN2拷贝数，预期的scale_factor应接近0.75（2:1）或0.5（2:0），在这两种情况下，每个位置的比例Pi_p应接近1/2。
#===============================================================================================================
Pi_a=data["Pi_a"][sample]
Pi_b=data["Pi_b"][sample]
Pi_c=data["Pi_c"][sample]
carrier_cutoff_pi=0.33333333
carrier_silence_p1=data["g.27134T>G"][sample].split(" ")[0]
carrier_silence_p2=data["g.27706_27707delAT"][sample].split(" ")[0]
bing=0.1
#3.1 silent carriers pi~=0.5:muts and check scale factor ~=0.75 or 0.5.
if (carrier_silence_p1 != "T") and (carrier_silence_p2 != "AT"):
    workout="SMN1_2+0_pred"
#3.2 normal pi>=0.3
elif (Pi_a >= carrier_cutoff_pi) and (Pi_b >= carrier_cutoff_pi) and (Pi_c >= carrier_cutoff_pi):
    workout="SMN1_exon7_normal|纯合|2/2"
#3.3 carrier pi< 0.3
elif (bing < Pi_a < carrier_cutoff_pi) or (bing < Pi_b < carrier_cutoff_pi) or (bing < Pi_c < carrier_cutoff_pi):
    workout="SMN1_exon7_del|杂合|2/1" 
#3.4 bing pi<=0.1
elif (Pi_a <= bing) or (Pi_b <= bing) or (Pi_c <= bing):
    workout="SMN1_exon7_del|纯合|2/0" 
#3.5 other bug
else:
    print("SMA workout error.")
with open(sma_workout_outputf_wd, 'w') as file:
    file.write(str(workout) + '\n')
