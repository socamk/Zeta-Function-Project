from mpmath import iv, isnan
import argparse
import sys
from old_verification.interval_arithmetic import my_atan

def read_zeros(file_name):
    '''
    This function takes a file containing zeros of the Riemann zeta function
    and adds them to a list

    input: the name of a text file containing zeros of the zeta function
    
    output: list containing those zeros
    '''
    file = open(file_name)      #open the text file
    zeros = []                 #create an empty list to hold the zeros
    iv.dps = 15
    for line in file:           #loop through the lines of the file
        words = line.split()    #split the line on whitespace and add the strings to a list
        zero = iv.mpf(words[1])  #turn the string into an interval    
        zeros.append(zero)   #add the interval to the list
    return zeros  


def cl_zero(point, power, gamma):
    match point:
        case 1:
            if power == 1:
                return iv.mpf("1")/ (iv.mpf("0.25") + gamma**iv.mpf("2"))
            elif (power % 1) == 0:
                rho = complex(iv.mpf(0.5), gamma)
                total = iv.mpf(1) / (rho ** power)
                return iv.mpf(2) * total.real
            else:
                sys.exit("Invalid power")
        case 0.5:
            if (power % 2 == 0):
                numerator = iv.mpf("2") * (iv.mpf("-1")**(power/iv.mpf("2")))
                denominator = gamma ** (power)
                return numerator/denominator
            else:
                sys.exit("Invalid power")
        case _:
            sys.exit("Invalid expansion point")

def ce_zero(point, power, gamma):
    match point:
        case 1:
            if power == 1:
                return iv.mpf("2")/ (iv.mpf("1") + gamma**iv.mpf("2"))
            else:
                sys.exit("Invalid power")
        case 0.5:
            if (power % 2 == 0):
                numerator = iv.mpf("2") * gamma
                numerator = my_atan(numerator)
                numerator = iv.mpf(power) * numerator
                numerator = iv.cos(numerator)
                numerator = iv.mpf([iv.absmin(numerator), iv.absmax(numerator)])
                numerator = numerator * iv.mpf("4")
                denominator = (iv.mpf("0.25") + (gamma ** iv.mpf("2"))) ** (iv.mpf(power)/iv.mpf("2"))
                return numerator / denominator
            else:
                sys.exit("Invalid power")

def general_sum(point, power, zeros):
    sum = iv.mpf("0")
    for zero in zeros:
        term = cl_zero(point, power, zero)
        sum += term
    return sum 


def verify_RH_list(zeros, heights, maximum, tail, power, point):
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
    t0_contribution = ce_zero(point, power, t0)       #find the amount that t0 adds to the sum
    sum = iv.mpf("0")         #initialize the sum
    results = []        #initialize an empty list to hold the results
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = general_sum(point, power, zeros)
    for zero in zeros:      #loop through all the zeros given
        next_term = cl_zero(point, power, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        if (sum + t0_contribution) > iv.mpf(maximum):       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            t0 = next(heights)              #set t0 to the next height in the list
            t0_contribution = ce_zero(point, power, t0)       #find contribution of the next height
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
    t0 = iv.mpf(start)      #set t0 to the given initial value
    t0_contribution = value(iv.mpf("1.0"), t0)       #find the amount that t0 adds to the sum
    sum = iv.mpf("0")         #initialize the sum
    results = []        #initialize an empty list to hold the results
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = full_sum(zeros)
    for zero in zeros:      #loop through all the zeros given
        next_term = value(iv.mpf("0.5"), zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        while (sum + (t0_contribution)) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            t0 += iv.mpf(step)      #increment t0 by the amount given
            t0_contribution = value(iv.mpf("1.0"), t0)       #find contribution of the new height
        index += 1      #increment the index at the end of each loop
    return results      #when entire loop is finished, return the results


def main():
    #make ArgumentParser object
    parser = argparse.ArgumentParser(description='Count the number of zeros needed to verify the Riemann Hypothesis underneath a series of heights.')
    #create all of the options that the program can use
    parser.add_argument('-i', '--interval', action='store', nargs=2, type=float, help='interval of heights to use when verifying the RH', metavar=('START', 'STEP'))
    parser.add_argument('-z', '--zeros', action='store', type=int, help='maximum number of zeros to use')
    parser.add_argument('-m', '--method', action='store', default=[1, 1], help='method to use for verification. First argument is the expansion point and second is the power to use.', type=float, nargs=2)
    parser.add_argument('-p', '--precision', action='store', nargs='+', help='precisions to use for the infinite sum', default=[2,3], type=int)
    parser.add_argument('-t', '--tail', action='store_true', help='include the sum of the tail in the results')
    parser.add_argument('-s', '--sum', action='store_true', help='take the sum over the zeros instead of verifying the RH')
    parser.add_argument('-f', '--file', action='store', type=str, help='alternate file of zeros to use')
    #process all of the options given in the command line arguments
    args = parser.parse_args()
    iv.dps = 20
    iv.pretty = True
    if args.file == None:
        all_zeros = read_zeros("zeros.txt")
    else:
        try:
            all_zeros = read_zeros(args.file)
        except:
            sys.exit("Could not read file, please try again")    
    if args.zeros == None:
        zeros = all_zeros
    else:
        zeros = all_zeros[:args.zeros]
    #determine method to use in the verifcation function (only one option right now)
    match args.method[0]:
        case 1:
            method_string = "1/rho"
            all_sums = ["0.02309570896612103381"]
        case 0.5:
            method_string = "1/(0.5 - rho)^" + str(args.method[1])
            all_sums = ["NaN","NaN","NaN","0.000074345198621964", "NaN", "-0.000000288347862801946559390523"]
        case _:
            sys.exit("Invalid expansion point")
    #the values of different methods to different precisions which can be used in the verification function
    #only one sublist because only one method exists at the moment
    #3 and 4 significant figure values are the same because they round tp the same amount
    if args.sum == True:
        for prec in args.precision:
            iv.dps = prec - 5
            new_zeros = []
            for zero in zeros:
                new_zeros.append(iv.mpf(zero))
            print(general_sum(args.method[0], args.method[1], new_zeros))
    else:
        #determine which sums are going to be used based on method and precision given
        sums = []
        for choice in args.precision:
            #iv.dps = choice - 5
            try:
                num = iv.mpf(all_sums[int(args.method[1] - 1)])
                sums.append(num)
                if isnan(num.a):
                    raise NotImplementedError
            except:
                print(all_sums)
                print(args.precision)
                print(args.method)
                sys.exit("Sorry, this power has not been implemented yet")
        #for each sum, verify using the interval if given, the list of zeros if not
        iv.dps = 30
        for sum in sums:
            if args.interval == None:
                print("here")
                results = verify_RH_list(zeros, zeros, sum, args.tail, args.method[1], args.method[0])
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
                    print(number, end="\t")
                print()
            print()

if __name__ == '__main__':
    main()