import argparse

def read_zeros(file_name):
    '''
    This function takes a file containing zeros of the Riemann zeta function
    and adds them to a list

    input: the name of a text file containing zeros of the zeta function
    
    output: list containing those zeros
    '''
    file = open(file_name)      #open the text file
    zeros = []                 #create an empty list to hold the zeros
    for line in file:           #loop through the lines of the file
        words = line.split()    #split the line on whitespace and add the strings to a list      
        zeros.append(float(words[1]))   #take the second element of the list, cast it to a float, and add it to the list
    return zeros        



def inverse_rho(beta, gamma):
    '''
    This function calculates 1/rho for a given beta and gamma, 
    where it is assumed that rho = beta + gamma*i

    input: 
        beta: the real part of a zero of the zeta function
        gamma: the imaginary part of a zero of the zeta function

    output: 1/(beta + gamma*i)
    '''
    return (beta/((beta**2) + (gamma**2)))


def full_sum(zeros, value):
    '''
    This function finds the total sum of terms given a list of zeros and
    a method to use on each term

    input:
        zeros: list of imaginary parts of the Riemann zeros
        value: function to use for verification (1/rho, 1/rho^2, etc)

    output:
        the sum of all terms
    '''
    sum = 0.0
    for zero in zeros:
        term = 2 * value(0.5, zero)
        sum = sum + term
    return sum
        


def verify_RH_list(zeros, heights, value, maximum, tail):
    '''
    This function finds the number of zeros needed to verify the Riemann Hypothesis
    at a series of heights. Currently, the only method it can use is the sum of 1/rho, 
    but it has flexibility, so that other powers of rho can be added when necessary.
    
    input: 
        zeros: list of imaginary parts of the Riemann zeros
        heights: list of the heights which the RH should be verified under
        value: function to use for verification (1/rho, 1/rho^2, etc)
        maximum: value of the sum that indicates a contradiction
        tail: whether to include the tail sum in the results

    output: list containing heights where the RH could be confirmed, the number of zeros required, and optionally the tail sum
    '''
    index = 0
    heights = iter(heights)     #make an iterator through the list of heights
    t0 = next(heights)      #set t0 to the first height in the list
    t0_contribution = value(1, t0)       #find the amount that t0 adds to the sum
    sum = 0         #initialize the sum
    results = []        #initialize an empty list to hold the results
    if tail:
        total_sum = full_sum(zeros, value)
    for zero in zeros:      #loop through all the zeros given
        next_term = 2 * value(0.5, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        if (sum + t0_contribution*2) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, add height and number of zeros to the results
            if tail:
                tail_sum = total_sum - sum
                result.append(tail_sum)
            results.append(result)
            t0 = next(heights)              #set t0 to the next height in the list
            t0_contribution = value(1.0, t0)       #find contribution of the next height
        index += 1      #increment the index at the end of each loop
    return results      #when entire loop is finished, return the results


def verify_RH_interval(zeros, start, step, value, maximum):
    '''
    Same as above, but this function evaluates at regular intervals instead of a list.
    
    input: 
        zeros: list of imaginary parts of the Riemann zeros
        start: first value to evaluate at
        step: length of interval to use
        value: function to use for verification (1/rho, 1/rho^2, etc)
        maximum: value of the sum that indicates a contradiction

    output: list containing heights where the RH could be confirmed and the number of zeros required
    '''
    index = 0
    t0 = start      #set t0 to the given initial value
    t0_contribution = value(1.0, t0)       #find the amount that t0 adds to the sum
    sum = 0         #initialize the sum
    results = []        #initialize an empty list to hold the results
    for zero in zeros:      #loop through all the zeros given
        next_term = 2 * value(0.5, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        if (sum + t0_contribution) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            results.append([t0, index + 1])     #if so, RH is verified, add height and number of zeros to the results
            t0 += step      #increment t0 by the amount given
            t0_contribution = value(1.0, t0)       #find contribution of the new height
        index += 1      #increment the index at the end of each loop
    return results      #when entire loop is finished, return the results




def main():
    parser = argparse.ArgumentParser(description='Count the number of zeros needed to verify the Riemann Hypothesis underneath a series of heights.')
    parser.add_argument('-i', '--interval', action='store', nargs=2, type=float, help='interval of heights to use when verifying the RH', metavar=('START', 'STEP'))
    parser.add_argument('-z', '--zeros', action='store', type=int, help='maximum number of zeros to use')
    parser.add_argument('-m', '--method', action='store', choices=[1], default=1, help='method to use for verification. Enter 1 for 1/rho, 2 for 1/rho^2, etc.', type=int)
    parser.add_argument('-p', '--precision', action='store', choices=[2, 3, 4], nargs='+', help='precisions to use for the infinite sum', default=[2,3], type=int)
    parser.add_argument('-t', '--tail', action='store_true', help='include the sum of the tail in the results')
    args = parser.parse_args()
    #zeros.txt contains the imaginary part of first 100,000 zeros of the Riemann zeta function
    #zeros.txt was sourced from lmfdb.org and information about these zeros can be found at
    #lmfdb.org/zeros/zeta
    all_zeros = read_zeros("zeros.txt")
    if args.zeros == None:
        zeros = all_zeros
    else:
        zeros = all_zeros[:args.zeros]
    if args.method == 1:
        method = inverse_rho
        method_string = "1/rho"
    all_sums = [[0.024, 0.0231, 0.02310]]
    sums = []
    for choice in args.precision:
        sums.append(all_sums[args.method - 1][choice - 2])
    for sum in sums:
        if args.interval == None:
            results = verify_RH_list(zeros, zeros, method, sum, args.tail)
        else:
            results = verify_RH_interval(zeros, args.interval[0], args.interval[1], method, sum)
        print("These results were generated using the sum of", method_string, "must be less than", sum)
        print("{:<20} {:<30}".format("Height", "Zeros needed to verify RH"), end='')
        if args.tail:
            print("Tail Sum")
        else:
            print()
        for result in results:
          print("{:<20} {:<30}".format(result[0], result[1]), end='')
          if args.tail:
            print(result[2])
        else:
            print()



if __name__ == '__main__':
    main()
