import random
import csv
import statistics
import array
#CAP280
#Pre-ART Decay, half life of 191 weeks
#Chi-Squared statistic
#expected seqs rounded to nearest whole number
#midpoint bins


def numbers():
    bin1 = 0
    bin2 = 0
    bin3 = 0
    bin4 = 0
    bin5 = 0
    bin6 = 0
    bin7 = 0
    bin8 = 0
    bin9 = 0
    sample = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    output = random.choices(sample, weights=(6, 7, 8, 9, 13, 15, 13, 10, 13, 5), k=39)
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
        if i == 6:
            bin6 = bin6 + 1
        if i == 7: 
            bin7 = bin7 + 1
        if i == 8: 
            bin8 = bin8 + 1
        if i == 9:
            bin9 = bin9 + 1
    a = [bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9] # the number of DNA sequences in each pre-ART bin after each iteration
    mean = statistics.mean(a)
    weights = [0.06, 0.07, 0.08, 0.09, 0.13, 0.15, 0.13, 0.1, 0.13, 0.05]
    early_seqs = 39
    expected_seqs = []
    obs_minus_expected = []
    for i in weights:
        expected_seqs.append((round(i * early_seqs))) #calculate the expected number of sequences in each bin, given the probabilities
    for i in range(0, len(a)): #calculate the (observed-expected)^2/expected value for each bin
        obs_minus_expected.append((a[i]-expected_seqs[i])**2/expected_seqs[i])
    sum_obs_minus_expected = sum(obs_minus_expected) #iteration sum of the difference between observed and expected sequence counts. 

    
    return(a, expected_seqs, obs_minus_expected, sum_obs_minus_expected)

print(numbers())

with open('240822_CAP280_ChiSquare_PreART_Decay_191week_halflife_Rounded_Midpoint_Bins.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['DNA Sequences per pre-ART Timepoint', 'Expected_Seqs', 'ChiSquared_Observed_Minus_Expected',\
             'Iteration_Sum_ChiSquared_O_E_Difference']
    writer.writerow(header)
    n = 0
    while n < 10000:
        writer.writerow(numbers())
        n = n + 1