import argparse
from decimal import *
import sys
import math
import cmath

def read_zeros(file_name):
    '''
    This function takes a file containing zeros of the Riemann zeta function
    and adds them to a list

    input: the name of a text file containing zeros of the zeta function
    
    output: list containing those zeros
    '''
    file = open(file_name)      #open the text file
    context = getcontext()
    context.prec = 100
    context.rounding = ROUND_UP
    zeros = []                 #create an empty list to hold the zeros
    for line in file:           #loop through the lines of the file
        words = line.split()    #split the line on whitespace and add the strings to a list      
        decimal = context.create_decimal(words[1])     #turn the string into a decimal
        zeros.append(decimal)   #add the decimal to the list
    return zeros        

def cl_zero_at_half(power, gamma):
    if (power % 2 == 0):
        numerator = Decimal("2") * (Decimal("-1")**(power/Decimal("2")))
        denominator = gamma ** (power)
        return numerator/denominator
    else:
        sys.exit("Invalid power")

def cl_zero_at_one(power, gamma):
    if power == 1:
        return Decimal("1")/ (Decimal("0.25") + gamma**Decimal("2"))
    elif (power % 1) == 0:
        rho = complex(Decimal(0.5), gamma)
        total = 1 / (rho ** power)
        return 2 * total.real
    else:
        sys.exit("Invalid power")

def ce_zero_at_half(power, gamma):
    if (power % 2 == 0):
        #numerator = abs(math.cos(Decimal("2") * power * math.atan(Decimal("2") * gamma))) * Decimal("4")
        numerator = Decimal("2") * gamma
        numerator = math.atan(numerator)
        numerator = Decimal(numerator)
        numerator = power * numerator
        numerator = math.cos(numerator)
        numerator = Decimal(numerator)
        numerator = abs(numerator)
        numerator = numerator * Decimal("4")
        
        denominator = (Decimal("0.25") + (gamma ** Decimal("2"))) ** (power/Decimal("2"))
        return numerator/denominator
    else:
        sys.exit("Invalid power")

def ce_zero_at_one(power, gamma):
    if power == 1:
        return Decimal("2")/ (Decimal("1") + gamma**Decimal("2"))
    else:
        sys.exit("Invalid power")

def general_sum(zeros, power, cl_method):
    sum = Decimal("0")
    for zero in zeros:
        term = cl_method(power, zero)
        sum += term
    return sum

        


def verify_RH_list(zeros, heights, maximum, tail, cl_method, ce_method, power):
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
    power = Decimal(power)
    heights = iter(heights)     #make an iterator through the list of heights
    t0 = next(heights)      #set t0 to the first height in the list
    t0_contribution = ce_method(power, t0)       #find the amount that t0 adds to the sum
    sum = Decimal("0")         #initialize the sum
    results = []        #initialize an empty list to hold the results
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = general_sum(zeros, power, cl_method)
    for zero in zeros:      #loop through all the zeros given
        next_term = cl_method(power, zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        while (abs(sum) + abs(t0_contribution)) >= abs(maximum):       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            t0 = next(heights)              #set t0 to the next height in the list
            t0_contribution = ce_method(power, t0)       #find contribution of the next height
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
    start = str(start)
    t0 = Decimal(start)      #set t0 to the given initial value
    t0_contribution = value(Decimal("1.0"), t0)       #find the amount that t0 adds to the sum
    sum = Decimal("0")         #initialize the sum
    results = []        #initialize an empty list to hold the results
    if tail:        #if including tail sum values, find the value of the sum over all zeros given
        total_sum = full_sum(zeros, value)
    for zero in zeros:      #loop through all the zeros given
        next_term = value(Decimal("0.5"), zero)        #find how much the current zero will add to the sum
        sum = sum + next_term           #add the amount to the sum
        while (sum + t0_contribution) > maximum:       #see if the sum plus the contribution from t0 exceeds the maximum
            result = [t0, index + 1]     #if so, RH is verified, record height and number of zeros
            if tail:        #if including tail sum values, find the value by subtracting the current sum from the total sum
                tail_sum = total_sum - sum
                result.append(tail_sum)     #add the tail sum to the result list
            results.append(result)      #add all the results to the results list
            step = str(step)
            t0 += Decimal(step)      #increment t0 by the amount given
            t0_contribution = value(Decimal("1.0"), t0)       #find contribution of the new height
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
    #zeros.txt contains the imaginary part of first 100,000 zeros of the Riemann zeta function
    #zeros.txt was sourced from lmfdb.org and information about these zeros can be found at
    #lmfdb.org/zeros/zeta
    if args.file == None:
        all_zeros = read_zeros("zeros.txt")
    else:
        try:
            all_zeros = read_zeros(args.file)
        except:
            sys.exit("Could not read file, please try again")    
    #determine number of zeros to use in the verifcation function
    if args.zeros == None:
        zeros = all_zeros
    else:
        zeros = all_zeros[:args.zeros]
    #determine method to use in the verifcation function (only one option right now)
    getcontext().rounding = ROUND_UP
    match args.method[0]:
        case 1:
            cl_method = cl_zero_at_one
            ce_method = ce_zero_at_one
            method_string = "1/rho"
            #the values of different methods to different precisions which can be used in the verification function
            all_sums = ["0.02309570896612103381"]
        case 0.5:
            cl_method = cl_zero_at_half
            ce_method = ce_zero_at_half
            method_string = "1/(0.5 - rho)^" + str(args.method[1])
            #S: -0.0000002883482314549231
            all_sums = ["NaN","NaN","NaN","0.000074345198621964", "NaN", "-0.000000288347862801946559390523"]
            sums = []
        case _:
            sys.exit("Invalid expansion point")
    if args.sum == True:
        for prec in args.precision:
            getcontext().prec = prec
            getcontext().rounding = ROUND_UP
            new_zeros = []
            for zero in zeros:
                new_zeros.append(getcontext().create_decimal(zero))
            getcontext().rounding = ROUND_DOWN
            print(general_sum(new_zeros, Decimal(args.method[1]), cl_method))
    else:
        #determine which sums are going to be used based on method and precision given
        sums = []
        for choice in args.precision:
            getcontext().prec = choice
            try:
                num = getcontext().create_decimal(all_sums[int(args.method[1] - 1)])
                sums.append(num)
                if num.is_nan():
                    raise NotImplementedError
            except:
                sys.exit("Sorry, this power has not been implemented yet")
        #for each sum, verify using the interval if given, the list of zeros if not
        getcontext().rounding = ROUND_DOWN
        getcontext().prec = 100
        for sum in sums:
            if args.interval == None:
                results = verify_RH_list(zeros, zeros, sum, args.tail, cl_method, ce_method, args.method[1])
            else:
                results = verify_RH_interval(zeros, args.interval[0], args.interval[1], method, Decimal(sum), args.tail)
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
                    getcontext().prec = 15
                    number = getcontext().create_decimal(number)
                    print("{:<30}".format(number), end="")
                print()
            print()



if __name__ == '__main__':
    main()
