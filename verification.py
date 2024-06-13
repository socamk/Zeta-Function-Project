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


def inverse_fourth_power(beta, gamma):
    return (-6.0/(gamma**4))


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
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = full_sum(zeros, value)
    for zero in zeros:      #loop through all the zeros given
        next_term = 2 * value(0.5, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        while (sum + t0_contribution*2) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            t0 = next(heights)              #set t0 to the next height in the list
            t0_contribution = value(1.0, t0)       #find contribution of the next height
        index += 1      #increment the index at the end of each loop
    return results      #when entire loop is finished, return the results


def verify_RH_interval(zeros, start, step, value, maximum, tail):
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
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = full_sum(zeros, value)
    for zero in zeros:      #loop through all the zeros given
        next_term = 2 * value(0.5, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        while (sum + t0_contribution) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            t0 += step      #increment t0 by the amount given
            t0_contribution = value(1.0, t0)       #find contribution of the new height
        index += 1      #increment the index at the end of each loop
    return results      #when entire loop is finished, return the results




def main():
    #make ArgumentParser object
    parser = argparse.ArgumentParser(description='Count the number of zeros needed to verify the Riemann Hypothesis underneath a series of heights.')
    #create all of the options that the program can use
    parser.add_argument('-i', '--interval', action='store', nargs=2, type=float, help='interval of heights to use when verifying the RH', metavar=('START', 'STEP'))
    parser.add_argument('-z', '--zeros', action='store', type=int, help='maximum number of zeros to use')
    parser.add_argument('-m', '--method', action='store', choices=[1], default=1, help='method to use for verification. Enter 1 for 1/rho, 2 for 1/rho^2, etc.', type=int)
    parser.add_argument('-p', '--precision', action='store', choices=[2, 3, 4, 5, 6, 7, 8], nargs='+', help='precisions to use for the infinite sum', default=[2,3], type=int)
    parser.add_argument('-t', '--tail', action='store_true', help='include the sum of the tail in the results')
    #process all of the options given in the command line arguments
    args = parser.parse_args()
    #zeros.txt contains the imaginary part of first 100,000 zeros of the Riemann zeta function
    #zeros.txt was sourced from lmfdb.org and information about these zeros can be found at
    #lmfdb.org/zeros/zeta
    all_zeros = read_zeros("zeros.txt")
    #determine number of zeros to use in the verifcation function
    if args.zeros == None:
        zeros = all_zeros
    else:
        zeros = all_zeros[:args.zeros]
    #determine method to use in the verifcation function (only one option right now)
    if args.method == 1:
        method = inverse_rho
        method_string = "1/rho"
    #the values of different methods to different precisions which can be used in the verification function
    #only one sublist because only one method exists at the moment
    #3 and 4 significant figure values are the same because they round tp the same amount
    all_sums = [[0.024, 0.0231, 0.02310, 0.023096, 0.0230958, 0.02309571, 0.023095709]]
    #determine which sums are going to be used based on method and precision given
    sums = []
    for choice in args.precision:
        sums.append(all_sums[args.method - 1][choice - 2])
    #for each sum, verify using the interval if given, the list of zeros if not
    for sum in sums:
        if args.interval == None:
            results = verify_RH_list(zeros, zeros, method, sum, args.tail)
        else:
            results = verify_RH_interval(zeros, args.interval[0], args.interval[1], method, sum, args.tail)
        #print the results of the verification
        print("These results were generated using the sum of", method_string, "must be less than", sum)
        result_strings = ["Height", "Zeros needed to verify RH"]
        if args.tail:
            result_strings.append("Tail Sum")
        for word in result_strings:
            print("{:<30}".format(word), end="")
        print()
        for result in results:
            for number in result:
                print("{:<30}".format(number), end="")
            print()
        print()



if __name__ == '__main__':
    main()
