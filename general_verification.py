import subprocess, bisect, sys, argparse
from mpmath import iv, nprint, nstr
from enum import Enum
from tail_approximation import r, R


class Function(Enum):
    '''
    Enum to keep track of the different L-functions that can be verified
    '''
    RIEMANN = 1
    REAL_DIRICHLET = 2
    RAMANUJAN = 3
    ELLIPTIC = 4

class Verification(Enum):
    '''
    Enum to keep track of the different types of verification
    '''
    RIEMANN_HYPOTHESIS = 1
    COMPLETENESS = 2



def read_hiary_zeros(start, file_name, lines):
    '''
    Function to read in zeros using the format found on Dr. Ghaith Hiary's webpage

        start: value to add to each zero, should be a power of 10
        file_name: name of a file containing imaginary parts of the zeros
        lines: number of lines to read, creating intervals can take a while, so this option
            allows the user to only read in a subset of lines without modifying the file
    '''
    file = open(file_name)      #open the text file
    zeros = []                 #create an empty list to hold the zeros
    start = iv.mpf(start)       #create an interval with the starting value
    error = iv.mpf("1e-10")     #create an interval with the error bounds, taken from Dr. Hiary's webpage
    for j in range(lines):           #for loop determines how many lines will be read
        line = file.readline()      #read a line from the file
        words = line.split()    #split the line on whitespace and add the strings to a list
        if len(words) > 1:      #ensures we have not reached the end of the file
            word = words[0] + words[1][1:]      #concatenate the strings
            num = iv.mpf(word)      #turn the string into an interval
            zero = iv.mpf([num.a - error, num.b + error])   #adjust the interval by the error bounds
            zero = zero + start     #add the starting value to the interval
            zeros.append(zero)      #add the interval to the list of zeros
    return zeros  

def read_zeros(file_name, index):
    '''
    This function takes a file containing zeros of the Riemann zeta function
    and adds them to a list

    input: the name of a text file containing zeros of the zeta function
    
    output: list containing those zeros
    '''
    file = open(file_name)      #open the text file
    zeros = []                 #create an empty list to hold the zeros
    error = iv.mpf("1e-8")
    for line in file:           #loop through the lines of the file
        words = line.split()    #split the line on whitespace and add the strings to a list
        zero = iv.mpf(words[index])  #turn the string into an interval    
        if zero != iv.mpf("0"):
            zero = iv.mpf([zero.a - error, zero.b + error])
        zeros.append(zero)   #add the interval to the list
    return zeros  

def find_closest_index(zeros, y):
    '''
    Internal function to find the starting index in a list of zeros

    inputs: 
        zeros - list of imaginary parts of zeros of an L-Function
        y - imaginary part of the expansion point, should be input as a string

    output: list containing index and value of the largest zero with imaginary part less than the expansion point
    '''
    y = iv.mpf(y)       #make sure y is an interval
    length = len(zeros)     #get the length of the list of zeros
    index = bisect.bisect_left(zeros, y)    #find index where y would be inserted, see python documentation for more information
    #check edge cases
    if index == 0:
        if y in zeros[0] and y > iv.mpf("0"):
            sys.exit("Invalid expansion point, please try again")
        return [0, zeros[0]]
    elif index == length:
        return [length - 1, zeros[length - 1]]
    else:
        #get values to the left and right of y
        left_val = zeros[index - 1]
        right_val = zeros[index]
        #ensure that y is not contained in either of those values
        if y in left_val or y in right_val:
            sys.exit("Invalid expansion point, please try again")
        #return index and value that are left of y
        return [index - 1, left_val]


def von_mangoldt_term(N, x, y, function, input):
    '''
    Internal function to calculate the part of the sum involving the Von Mangoldt function
    inputs:
        x - string, real part of the expansion point
        y - string, imaginary part of the expansion point
        function - enum for the type of function being evaluated
        file name - name of a file containing e^Lambda(n) or Lambda(n)
    
    output: interval containing this portion of the sum using N terms for the sum involving the
    Von Mangoldt function

        TODO - handle case where N is greater than available VM terms
    '''
    if (isinstance(input, str)):
        file = open(input)  #open file
        sum = iv.mpc("0")       #initiate sum
        if function.value == Function.RIEMANN.value: 
            for i in range(N):      #for loop determines how many terms will be used
                line = file.readline()      #read a line from the file
                if line.strip() != "1":     #if line does not equal 1, meaning log(line) != 0
                    sum += (iv.log(iv.mpf(line.strip())) / (iv.mpf(i + 1) ** (iv.mpc("1", "0") - iv.mpc(x, y))))    #use the line to calculate the next term and add it to the sum
            sum = (iv.mpc("-1","0") * sum) - (iv.mpc("1","0") / iv.mpc(x, y)) #multiply sum by -1 and subtract 1/z
        elif function.value >= Function.RIEMANN.value:    #for functions other than zeta
            #ensure a real expansion point is being used
            if y != "0":
                sys.exit("Invalid expansion point. Please choose a point on the real line and try again.")
            #same as Riemann case, loop through the file
            for i in range(1, N + 1):
                line = file.readline()  #read a line from the file
                #if the line is not zero, calculate the next term and add it to the sum
                if line.strip() != "0":
                    sum += iv.mpf(line.strip())/(iv.mpf(i) ** (iv.mpf("1") - iv.mpf(x)))
            #multiply the sum by -1
            sum = iv.mpf("-1") * sum
        return sum
    elif (isinstance(input, list)):
        value = iv.mpf(input[0])
        error = iv.mpf(input[1])
        sum = iv.mpf([value.a - error.b, value.b + error.b])
        return sum
    else:
        sys.exit("Invalid input, please try again.")


def error_term(N, x, function):
    '''
    Internal function to calculate the truncation error from the sum over primes

    input:
        N - number of terms used in the sum over primes
        x - real part of the expansion point
        function - enum representing the type of L-function being evaluated
    '''
    #set the value of r depending on what type of function is being evaluated
    if function.value > 2: 
        r = 2
    else:
        r = 1
    #ensure that N and x are intervals
    N = iv.mpf(N)
    x = iv.mpf(x)
    #calculate the error term, separated into three steps for readability
    first_term = (r * (N ** x)/x)
    second_term = iv.mpf("2.85") * (((iv.mpf("2") * x) - iv.mpf("1"))/iv.log(N))
    error = first_term * (second_term - iv.mpf("1"))
    return error

def digamma_term(x, y, function, d):
    '''
    TODO - clean this up and fix docs
    Internal function to calculate the portion of the sum involving the digamma function

    inputs:
        x - string, real part of the expansion point
        y - string, imaginary part of the expansion point
        function - enum representing the type of function being evaluated
        d - fundamental discriminant, input is None if not applicable

    output: interval containing (1/2)*Digamma(3/2 - z/2)

    '''
    #use command line to run compiled C program with two arguments and pipe stdout to a text file
    if function.value == Function.RIEMANN.value:
        #run program using FLINT for calculations
        process = subprocess.run(["./riemann_digamma", x, y], capture_output=True, encoding="utf-8") 
        #save output of the program
        line = process.stdout
        #split the line on whitespace and collect in a vector
        words = line.split("\n")
        values = []     #empty list to hold results
        #processing FLINT output into a string that mpmath can understand
        for word in words:
            nums = word.split(" +/- ")
            for i in range(len(nums)):
                nums[i] = nums[i].strip("[]")
            #if FLINT output has an error term, construct an interval using those error bounds
            if len(nums) == 2:
                base = iv.mpf(nums[0])
                error = iv.mpf(nums[1])
                upper = base.b + error.b
                lower = base.a - error.b
                values.append(iv.mpf([lower, upper])) #add to result list
            #if FLINT output has no error (usually if result is exactly zero), construct interval from the output value
            elif len(nums) == 1:
                values.append(iv.mpf(nums[0]))  #add to result list
        #result list should have two entries representing the real and imaginary parts of the calculation
        #create a complex interval in mpmath using those entries
        value = iv.mpc(values[0], values[1])
    elif function.value == Function.REAL_DIRICHLET.value:
        #set the value of m based on the sign of the conductor
        if iv.mpf(d) > iv.mpf("0"):
            m = "0"
        else:
            m = "1"
        #run program using FLINT for calculations
        process = subprocess.run(["./general_digamma", x, m], capture_output=True, encoding="utf-8") 
        #process output of FLINT program back into mpmath interval
        line = process.stdout
        word = line.strip()
        nums = word.split(" +/- ")
        for i in range(len(nums)):
            nums[i] = nums[i].strip("[]")
        base = iv.mpf(nums[0])
        error = iv.mpf(nums[1])
        upper = base.b + error.b
        lower = base.a - error.b
        value = iv.mpf([lower, upper])
    elif function.value >= Function.RAMANUJAN.value:
        #set values of m depending on the function
        if function.value == Function.RAMANUJAN.value:
            m1 = "5.5"
            m2 = "6.5"
        if function.value == Function.ELLIPTIC.value:
            m1 = "0.5"
            m2 = "1.5"
        value = iv.mpf("0")     #initialize the sum
        #for each m, find the digamma value using FLINT and add to the sum
        for val in [m1, m2]:
            process = subprocess.run(["./general_digamma", x, val], capture_output=True, encoding="utf-8")
            line = process.stdout 
            word = line.strip()
            nums = word.split(" +/- ")
            for i in range(len(nums)):
                nums[i] = nums[i].strip("[]")
            base = iv.mpf(nums[0])
            error = iv.mpf(nums[1])
            upper = base.b + error.b
            lower = base.a - error.b
            value += iv.mpf([lower, upper])
    #divide the final term by 2 and return it
    return iv.mpf("1/2") * value

def find_sum(x, y, N, function, d, file_name="Lambda_Values/Riemann_Lambda.txt"):
    '''
    Function to find the actual value of a sum over all zeros of the Riemann zeta function

    inputs:
        x - real part of the expansion point
        y - imaginary part of the expansion point
        N - number of terms to use to calculate the sum over Von Mangoldt values
        function - enum representing the type of function being evaluated
        d - fundamental discriminant, passed as None if not applicable
        file name - name of a file containing e^Lambda(n) for zeta and Lambda(n) for other functions

    output: interval containing the sum of 1/(rho - z) for all rho, using z = x + iy
    '''
    #find the logarithmic term of the sum based on the chosen function
    if function.value == Function.RIEMANN.value:
        log_term = iv.mpf("-1/2") * iv.log(iv.pi)
    elif function.value == Function.REAL_DIRICHLET.value:
        log_term = (iv.mpf("1/2") * iv.log(abs(int(d)))) - (iv.mpf("1/2") * iv.log(iv.pi))
    elif function.value == Function.RAMANUJAN.value:
        log_term = log_term = (iv.mpf("1/2") * iv.log(iv.mpf("1"))) - (iv.log(iv.pi))
    elif function.value == Function.ELLIPTIC.value:
        log_term = log_term = (iv.mpf("1/2") * iv.log(iv.mpf("37"))) - (iv.log(iv.pi))
    #find the term of the sum involving the digamma function
    dg_term = digamma_term(x, y, function, d)
    #find the term of the sum involving the sum over the primes
    vm_term = von_mangoldt_term(N, x, y, function, file_name)
    #find the truncation error from the sum over the primes
    e_term = error_term(N, x, function)
    #find the upper and lower bounds of the sum
    upper = log_term + dg_term + vm_term + e_term
    lower = log_term + dg_term + vm_term - e_term
    #return interval using those bounds
    return iv.mpf([lower.a, upper.b])


def ce_contribution(x, beta, eta):
    x = iv.mpf(x)
    eta = iv.mpf(eta)
    beta = iv.mpf(beta)
    num1 = beta - x
    den1 = (beta - x) ** 2
    den1 = den1 + (eta ** 2)
    val1 = num1/den1
    num2 = iv.mpf("1") - beta - x
    den2 = (iv.mpf("1") - beta - x) ** 2
    den2 = den2 + (eta ** 2)
    val2 = num2/den2
    return val1 + val2


def sum_over(zeros, x, y, function):
    '''
    Function to find the total contribution of a set of zeros of the Riemann Zeta Function

    input:
        zeros - set of intervals containing zeros of the zeta function
        x - real part of the expansion point
        y - imaginary part of the expansion point
        function - enum representing the function being evaluated
    output: interval representing the bounds of the sum contribution of the given zeros
    '''
    sum = iv.mpf("0")   #initialize sum
    x = iv.mpf(x)       #turn x and y into intervalz
    y = iv.mpf(y)
    beta = iv.mpf("1/2")    #set the interval for beta
    #if the function is zeta, each zero contributes (1/2 - x)/[(1/2 - x)^2 + (gamma - y)^2]
    if function.value == Function.RIEMANN.value:
        for zero in zeros:
            num = beta - x
            den = (beta - x) ** 2
            den = den + ((zero - y) ** 2)
            term = num/den
            sum += term
    #use different sum for a general L-function
    elif function.value >= Function.REAL_DIRICHLET.value:
        for zero in zeros:
            if zero == iv.mpf("0"):
                num = beta - x
                den = ((beta - x) ** iv.mpf("2"))
                sum += num/den
            else:
                num = iv.mpf("1") - (iv.mpf("2") * x)
                den = (iv.mpf("1/2") - x) ** iv.mpf("2")
                den = den + (zero ** iv.mpf("2"))
                term = num/den
                sum += term
    return sum

def verify(zeros, x, y, N, Tau, function, file, verification, tail=False, d=None):
    '''
    Function to verify a general L-function
    '''
    if float(x) >= 0:
        sys.exit("Innappropriate expansion point. Please choose a value of x < 0 and try again.")
    if float(y) < 0:
        sys.exit("Innappropriate expansion point. Please choose a value of y >= 0 and try again.")
    #find list of zeros inside the range given by tau
    upper_val = iv.mpf(y) + iv.mpf(Tau)
    lower_val = iv.mpf(y) - iv.mpf(Tau)
    first = find_closest_index(zeros, lower_val)
    last = find_closest_index(zeros, upper_val)
    if first[0] == 0 and first[1] > iv.mpf(y) - iv.mpf(Tau):
        zeros = zeros[:last[0] + 1]
    else:
        zeros = zeros[first[0] + 1:last[0] + 1]
    upper_bound = find_sum(x, y, N, function, d, file).real
    base_sum = sum_over(zeros, x, y, function)
    if verification == Verification.COMPLETENESS and tail == True:
        upper_tail_bound = R(x, y, Tau)
        total = base_sum + upper_tail_bound
        if total.b < upper_bound.a:
            print("The list given is incomplete")
            return
    #find bound on tail contribution if applicable
    if tail and function == Function.RIEMANN:
        lower_tail = r(x, y, Tau)
        base_sum = base_sum + lower_tail
    #check counterexamples and loop until no contradiction is reached
    done = False
    i = 1
    while not done:
        if verification == Verification.RIEMANN_HYPOTHESIS:
            val1 = ce_contribution(x, "1/2", i)
            val2 = ce_contribution(x, "1", i)
        elif verification == Verification.COMPLETENESS:
            val1 = ce_contribution(x, "1/2", i) * iv.mpf("1/2")
            val2 = ce_contribution(x, "0", i)
        contribution = min(val1, val2)
        if function.value >= Function.REAL_DIRICHLET.value:
            contribution = contribution * 2
        total = base_sum + contribution
        if total >= upper_bound:
            i += 1
        else:
            done = True
    #return largest integer that causes a contradiction
    return i - 1



def main():
    iv.dps = 40
    parser = argparse.ArgumentParser(description="Program to verify the Riemann Hypothesis or completeness within a subsection of a given list of zeros. Currently works with the Riemann zeta function, real Dirichlet functions, the Ramanujan tau function, and elliptic curves.")
    parser.add_argument("-R", "--Riemann", action='store', nargs=1, help='verify the Riemann zeta function around a point z = x + iy using zeros in a range of [y - τ, y + τ]', metavar="TAU")
    parser.add_argument("-D", "--Dirichlet", action='store', nargs=2, type=int, help='verify a real Dirichlet function around a point z = x + iy using zeros in a range of [y - τ, y + τ]', metavar=("CONDUCTOR", "TAU"))
    parser.add_argument("-T", "--Ramanujan", action='store', help='verify the Ramanujan tau function around a point z = x + iy using zeros in a range of [y - τ, y + τ]', metavar="TAU")
    parser.add_argument("-p", "--point", nargs=2, help="Point where the expansion is centered, default is -1", default=["-1", "0"], metavar=("REAL", "IMAGINARY"))
    parser.add_argument("-l", "--Lambda", nargs=2, help="File containing e^Λ(n) for zeta or Λ(n) for other functions and number of terms to use for the sum over primes", metavar=("FILENAME", "TERMS"))
    parser.add_argument("-H", "--H_zeros", nargs=3, help="use file of zero ordinates created by Dr. Ghaith Hiary", metavar=("FILENAME", "SHIFT", "LINES"))
    parser.add_argument("-z", "--zeros", nargs=2, help="use file of zero ordinates", metavar=("FILE_NAME", "COLUMN"))
    parser.add_argument("-t", "--tail", action='store_true', help='include upper and lower bounds on the tail of the sum in the verification. Currently only works for the Riemann zeta function')
    parser.add_argument("-c", "--completeness", action="store_true", help="verify completeness of a list of zeros instead of the Riemann Hypothesis")
    args = parser.parse_args()
    if args.Lambda == None:
        sys.exit("No Lambda values provided, please try again.")
    if args.zeros != None and args.H_zeros != None:
        sys.exit("Too many zero files. Please try again and provide one file with all zero ordinates.")
    elif args.zeros == None and args.H_zeros == None:
        sys.exit("No files with zero ordinates provided, please try again")
    elif args.zeros != None:
        zeros = read_zeros(args.zeros[0], args.zeros[1])
    elif args.H_zeros != None:
        zeros = read_hiary_zeros(args.H_zeros[1], args.H_zeros[0], int(args.H_zeros[2]))
    count = 0
    for arg in [args.Riemann, args.Ramanujan, args.Dirichlet]:
        if arg != None:
            count += 1
    if count == 0:
        sys.exit("No function provide, please try again.")
    elif count > 1:
        sys.exit("Multiple functions provided, please try again.")
    verification = Verification.RIEMANN_HYPOTHESIS
    if args.completeness == True:
        verification = Verification.COMPLETENESS
    if args.Riemann != None:
        val = verify(zeros, args.point[0], args.point[1], int(args.Lambda[1]), args.Riemann[0], Function.RIEMANN, args.Lambda[0], verification, args.tail)
    elif args.Dirichlet != None:
        val = verify(zeros, args.point[0], args.point[1], int(args.Lambda[1]), args.Dirichlet[1], Function.REAL_DIRICHLET, args.Lambda[0], verification, False, args.Dirichlet[0])
    elif args.Ramanujan != None:
        val = verify(zeros, args.point[0], args.point[1], int(args.Lambda[1]), args.Ramanujan[0], Function.RAMANUJAN, args.Lambda[0], verification)
    print("The list has been verified to a distance of", val)
if __name__ == "__main__":
    main()