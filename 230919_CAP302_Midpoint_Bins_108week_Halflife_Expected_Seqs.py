import random
import csv
import statistics
import array
#CAP302
#Pre-ART Decay, half life of 108 weeks
#Chi-Squared statistic
#expected seqs rounded to nearest whole number
#midpoint bins


def numbers():
    bin1 = 0
    bin2 = 0
    bin3 = 0
    bin4 = 0
    bin5 = 0
    sample = [1, 2, 3, 4, 5]
    output = random.choices(sample, weights=(10, 16, 23, 29, 22), k=32)
    for i in output:
        if i == 1:
            bin1 = bin1 + 1
        if i == 2:
            bin2 = bin2 + 1
        if i == 3:
            bin3 = bin3 + 1
        if i ==4:
            bin4 = bin4 + 1
        if i ==5:
            bin5 = bin5 + 1
    a = [bin1, bin2, bin3, bin4, bin5] # the number of DNA sequences in each pre-ART bin after each iteration
    mean = statistics.mean(a)
    weights = [0.10, 0.16, 0.23, 0.29, 0.22]
    early_seqs = 32
    expected_seqs = []
    obs_minus_expected = []
    for i in weights:
        expected_seqs.append((round(i * early_seqs))) #calculate the expected number of sequences in each bin, given the probabilities
    for i in range(0, len(a)): #calculate the (observed-expected)^2/expected value for each bin
        obs_minus_expected.append((a[i]-expected_seqs[i])**2/expected_seqs[i])
    sum_obs_minus_expected = sum(obs_minus_expected) #iteration sum of the difference between observed and expected sequence counts. 

    
    return(a, expected_seqs, obs_minus_expected, sum_obs_minus_expected)

print(numbers())

with open('230919_CAP302_ChiSquare_PreART_Decay_MidPoint_Bins_108week_halflife_Expected_Seqs_Rounded.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['DNA Sequences per pre-ART Timepoint', 'Expected_Seqs', 'ChiSquared_Observed_Minus_Expected',\
             'Iteration_Sum_ChiSquared_O_E_Difference']
    writer.writerow(header)
    n = 0
    while n < 10000:
        writer.writerow(numbers())
        n = n + 1
    