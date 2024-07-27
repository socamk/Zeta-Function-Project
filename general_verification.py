import subprocess, bisect, sys
from mpmath import iv, nprint, nstr
from enum import Enum
from tail_approximation import r, R


class Function(Enum):
    '''
    Enum to keep track of the different L-functions that can be verified
    '''
    RIEMANN = 1
    QUADRATIC = 2
    RAMANUJAN = 3
    ELLIPTIC = 4



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
    if function.value == 1:
        process = subprocess.run(["./riemann_digamma", x, y], capture_output=True, encoding="utf-8") 
    elif function.value == 2:
        if iv.mpf(d) > iv.mpf("0"):
            m = "0"
        else:
            m = "1"
        process = subprocess.run(["./general_digamma", x, m], capture_output=True, encoding="utf-8") 
    #get the output of the C program aptured by the Python subprocess
    if function.value == 1:
        line = process.stdout
        #split the line on whitespace and collect in a vector
        words = line.split("\n")
        values = []
        for word in words:
            nums = word.split(" +/- ")
            for i in range(len(nums)):
                nums[i] = nums[i].strip("[]")
            if len(nums) == 2:
                base = iv.mpf(nums[0])
                error = iv.mpf(nums[1])
                upper = base.b + error.b
                lower = base.a - error.b
                values.append(iv.mpf([lower, upper]))
            elif len(nums) == 1:
                values.append(iv.mpf(nums[0]))
        value = iv.mpc(values[0], values[1])
    elif function.value == 2:
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
    elif function.value >= 3:
        if function.value == 3:
            m1 = "5.5"
            m2 = "6.5"
        if function.value == 4:
            m1 = "0.5"
            m2 = "1.5"
        value = iv.mpf("0")
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
        file name - name of a file containing e^Lambda(n), defaults to existing file but 
        could also use a user-specified file

    output: interval containing the sum of 1/(rho - z) for all rho and z = x + iy
    '''
    if function.value == 1:
        log_term = iv.mpf("-1/2") * iv.log(iv.pi)
    elif function.value == 2:
        log_term = (iv.mpf("1/2") * iv.log(abs(int(d)))) - (iv.mpf("1/2") * iv.log(iv.pi))
    elif function.value == 3:
        log_term = log_term = (iv.mpf("1/2") * iv.log(iv.mpf("1"))) - (iv.log(iv.pi))
    elif function.value == 4:
        log_term = log_term = (iv.mpf("1/2") * iv.log(iv.mpf("37"))) - (iv.log(iv.pi))
    dg_term = digamma_term(x, y, function, d)
    print("dg =", dg_term)
    vm_term = von_mangoldt_term(N, x, y, function, file_name)
    print("vm =", vm_term)
    e_term = error_term(N, x, function)
    print("e =", e_term)
    print("log =", log_term)
    upper = log_term + dg_term + vm_term + e_term
    lower = log_term + dg_term + vm_term - e_term
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
        point - imaginary point to do the expansion at
    output: interval representing the bounds of the sum contribution of the given zeros
    '''
    sum = iv.mpf("0")
    x = iv.mpf(x)
    y = iv.mpf(y)
    beta = iv.mpf("1/2")
    if function.value == Function.RIEMANN.value:
        for zero in zeros:
            num = beta - x
            den = (beta - x) ** 2
            den = den + ((zero - y) ** 2)
            term = num/den
            sum += term
    elif function.value >= Function.QUADRATIC.value:
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

def verify(zeros, x, y, N, delta, function, file, type, tau=None, d=None):
    '''
    TODO - add method to include zeros by value instead of number
    TODO - make sure L-function inputs are correct
    '''
    if float(x) >= 0:
        sys.exit("Innappropriate expansion point. Please choose a value of x < 0 and try again.")
    if float(y) < 0:
        sys.exit("Innappropriate expansion point. Please choose a value of y >= 0 and try again.")
    mid = find_closest_index(zeros, y)
    print(mid)
    mid_index = mid[0]
    print(mid_index)
    length = len(zeros)
    if float(y) > 0:
        if delta > mid_index:
            sys.exit("Not enough zeros. Please use a smaller number of zeros or a higher expansion point and try again.")
        if delta + mid_index > length:
            sys.exit("Not enough zeros. Please use a smaller number of zeros or a lower expansion point and try again.")
        zeros = zeros[mid_index - delta: mid_index + delta]
    elif float(y) == 0:
        if delta + mid_index > length:
            sys.exit("Not enough zeros. Please use a smaller number of zeros or a lower expansion point and try again.")
        zeros = zeros[mid_index:mid_index + delta]
    else:
        sys.exit("Innappropriate expansion point. Please choose a value of y >= 0 and try again.")
    print(len(zeros))
    upper_bound = find_sum(x, y, N, function, d, file).real
    print("found sum =", upper_bound)
    base_sum = sum_over(zeros, x, y, function)
    print("sum over zeros =", base_sum)
    if tau != None:
        lower_tail = r(x, y, tau)
        print("tail =", lower_tail)
        base_sum = base_sum + lower_tail
    done = False
    i = 1
    last_val = None
    while not done:
        if type == 1:
            val1 = ce_contribution(x, "1/2", i)
            val2 = ce_contribution(x, "1", i)
        elif type == 2:
            val1 = ce_contribution(x, "1/2", i) * iv.mpf("1/2")
            val2 = ce_contribution(x, "0", i)
        contribution = min(val1, val2)
        if function.value >= Function.QUADRATIC.value:
            contribution = contribution * 2
        total = base_sum + contribution
        new_val = contribution
        if total >= upper_bound:
            #print(i, True)
            i += 1
            if i % 1000 == 0:
                print(i)
            last_val = new_val
        else:
            done = True
            #print("The Riemann Hypothesis has been verified between", float(y) - (i - 1), "and", float(y) + (i - 1))
            iv.dps = 10
            print("(" + x + "," + y + ")", end="\t")
            print(upper_bound, end="\t")
            print(base_sum, end="\t")
            print(last_val, end="\t")
            print(last_val + base_sum, end="\t")
            print(i - 1, end="\t")
            print()
            iv.dps = 40
    return i - 1



def main():
    iv.dps = 40

    
    #TEST CODE PLEASE IGNORE
    '''
    print(find_sum("-1", "0", 100000, Function.ELLIPTIC, None, ["0.547934606487896699144260368219", "1e-20"]))
    e_zeros = read_zeros("zeros/Elliptic_zeros.txt", 0)
    print(sum_over(e_zeros, "-1", "0", Function.ELLIPTIC))
    verify(e_zeros, "-1", "0", 100000, 100, Function.ELLIPTIC, ["0.547934606487896699144260368219", "1e-20"], 1)
    print(sum_over(e_zeros, "-1", "0", Function.ELLIPTIC) + (2 * ce_contribution("-1", "1", "10")))
    '''


    '''
    d_zeros = read_zeros("zeros/10^6_zeros.txt", 1)
    verify(d_zeros, "-1", "0", 100000, 5000, Function.QUADRATIC, "Lambda_Values/Lambda_Quadratic_-1159523.txt", 1, None, "-1159523")

    print(find_sum("-1", "0", 100000, Function.QUADRATIC, "-1159523", "Lambda_Values/Lambda_Quadratic_-1159523.txt"))
    print(sum_over(d_zeros, "-1", "0", Function.QUADRATIC))
    '''
    
    
    
    #zeros = read_hiary_zeros("1e28", "zeros/1e28.zeros.1000_10001001", 10000000)
    #y = iv.mpf("10000000000000000000000501675.8")
   # verify(zeros, "-2", "10000000000000000000000501675.8", 10000000, 4999999, Function.RIEMANN, "Lambda_Values/Riemann_Lambda.txt", 1, "501575.4")
    #verify(zeros, "-2", "10000000000000000000000501675.8", 10000000, 4999999, Function.RIEMANN, "Lambda_Values/Riemann_Lambda.txt", 2, "501575.4")

    '''
    print("R verification")
    r_zeros = read_zeros("zeros/Ramanujan_zeros.txt", 0)
    verify(r_zeros, "-1", "0", 100000, 20000, Function.RAMANUJAN, ["0.058326197419819564458801483366", "1e-20"], 1)
    print()
    print(find_sum("-1", "0", 100000, Function.RAMANUJAN, None, ["0.058326197419819564458801483366", "1e-20"]))
    sum = sum_over(r_zeros, "-1", "0", Function.RAMANUJAN)
    print("C =", sum)
    x = "-1"
    i = "84"
    val1 = 2 *min (ce_contribution(x, "1/2", i), ce_contribution(x, "1", i))
    print(sum + val1)
    i = "85"
    val2 = 2 * min (ce_contribution(x, "1/2", i), ce_contribution(x, "1", i))
    print(sum + val2)
    ''' 
    
    
    '''
    zeros = read_hiary_zeros("1e28", "zeros/1e28.zeros.1000_10001001", 10000000)
    
    total = iv.mpf(["31.4180626270347520984646803475482274464888561", "31.4180626270348458180158187325621374281859243"])
    sum = (sum_over(zeros, "-2", "10000000000000000000000501675.8", Function.RIEMANN))
    print(sum)
    tail = R("-2", "10000000000000000000000501675.8", "501575.4")
    print("R =", tail)
    print(sum + tail)
    print((sum + tail).b < total.a)
    print()
    val = zeros.pop(5200000)
    print("removed", val)
    sum = (sum_over(zeros, "-2", "10000000000000000000000501675.8", Function.RIEMANN))
    print("new sum =", sum)
    print(sum + tail)
    print((sum + tail).b < total.a)
    
    '''
    #verify(zeros, "-2", "10000000000000000000000501675.8", 10000000, 4999999, Function.RIEMANN, "Lambda_Values/Riemann_Lambda.txt", 1, "501575.4")
    #verify(zeros, "-2", "10000000000000000000000501675.8", 10000000, 4999999, Function.RIEMANN, "Lambda_Values/Riemann_Lambda.txt", 2, "501575.4")
    
if __name__ == "__main__":
    main()