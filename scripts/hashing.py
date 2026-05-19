import math


def Fj(j, m, d, epsilon):
    n = m * d / (1 + epsilon)
    if j <= 0 or j > m:
        return 0.0
    # Use log-gamma for log binomial coefficient to handle large values
    log_binom = math.lgamma(m + 1) - math.lgamma(j + 1) - math.lgamma(m - j + 1)

    
    # Compute ratio1 and its log term
    ratio1 = n * ((j / m) ** 2) / (j * d)
    if ratio1 <= 0:
        # return 0.0
        print("Ratio 1 less than 0")
    log_term2 = (j * d) * math.log(ratio1)


    # print("Hello")
    
    # Compute ratio2 and its log term

    if(n - j * d == 0): log_term3 = 0
    else:   
        ratio2 = (n - n * ((j / m) ** 2)) / (n - j * d)
        if ratio2 <= 0:  # Note: ratio2 should be <=1 in context, but guard
            # return 0.0
            print("Ratio 2 less than 0", n, j, m, d)
        log_term3 = (n - j * d) * math.log(ratio2)

    # Sum logs and exp
    log_Fj = log_binom + log_term2 + log_term3


    return math.exp(log_Fj)

def F(m, d, epsilon):
    n = m * d / (1 + epsilon)
    upper = math.floor(m / (1 + epsilon))
    total = 0.0
    for j in range(1, upper + 1):
        # print(j, m, d, n, epsilon)
        res = Fj(j, m, d, epsilon)
        total += res
    return total


def find_max_balls(m, d):

    res = 0
    for n in range(m * d - 10, d - 1, -1):
        epsilon = ((m * d) / n) - 1
        # print(epsilon)

        if(F(m, d, epsilon) < fail_prob): 
            res = n
            break

    return res


# Number of buckets
M = [8, 32, 128, 512, 2048]

# Capacity of each bucket
D = [8, 16]
# D = [4]


fail_prob = 2 ** (-40)


for m in M:
    for d in D:
        
        # For a given number of buckets m and bucket capacity d
        # find max balls n that can be allocated with fail_prob
        n = find_max_balls(m, d)
        print (m, d, n)


