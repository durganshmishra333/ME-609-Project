import random
import math
import matplotlib.pyplot as plt
import numpy as np


def plot_x_versus_iterations(x_array, iterations, ylabel, c, method_name):
    x = x_array
    k = iterations

    x_axis = []

    for i in range(0, k+1):
        x_axis.append(i)

    y_axis = x[abs(len(x_axis)-len(x)):]

    fig = plt.figure()
    axes = fig.add_subplot(111)

    plt.ylabel(ylabel)
    plt.xlabel("Number of iterations")

    axes.plot(x_axis, y_axis)

    for i in range(0, len(y_axis)):
        plt.plot(x_axis[i], y_axis[i], 'ro')

    plt.show(block=False)
    plt.savefig(f"range_plot_for_{ylabel}_{method_name}_part{c}.png")


def plot_range(a, b, new_a, new_b, c, d, method_name):
    x_axis = []
    y_axis = []

    for i in np.linspace(a, b, 1000):
        y = my_function(i, d, c)
        x_axis.append(i)
        y_axis.append(y)

    fig = plt.figure()
    axes = fig.add_subplot(111)

    plt.xlabel("x")
    plt.ylabel("f(x)")

    axes.plot(x_axis, y_axis)
    axes.plot([new_a], [my_function(new_a, d, c)], 'ro')
    axes.plot([new_b], [my_function(new_b, d, c)], 'ro')
    plt.show(block=False)
    plt.savefig(f"range_plot_{method_name}_part{c}.png")

# functions as given in the assignment


def my_function(x, d, c):
    if c == "1":
        eqn = (2*x-5)**4-(x**2-1)**3
    elif c == "2":
        eqn = 8 + x**3 - 2*x - 2 * math.exp(x)
    elif c == "3":
        eqn = 4 * x * math.sin(x)
    elif c == "4":
        eqn = 2 * (x-3)**2 + math.exp(0.5 * x**2)
    elif c == "5":
        eqn = x**2 - 10*math.exp(0.1*x)
    elif c == "6":
        eqn = 20*math.sin(x) - 15 * x**2
    if d == "N":
        return -1*eqn  # if the function needs to be maximized
    else:
        return eqn  # if the function needs to be minimized

# bounding phase function


def bounding_phase(a, b, d, c):
    org_a = a
    org_b = b
    r = random.uniform(a, b)  # assigning a random starting value
    # assigning a random value of delta
    delta = random.uniform(10**-9, 10**-12)

    K = 0
    # initial function evaluations
    F_P = my_function((r-delta), d, c)
    F_Q = my_function(r, d, c)
    F_R = my_function((r + delta), d, c)

    print("************")
    print("Bounding Phase Method")

    # reading out the data from iterations
    out = open(f"Bounding_Phase_iterations_part{c}.out", "w")

    # taking input till it satisfies the condition
    while 1:
        if F_P >= F_Q and F_Q >= F_R:
            delta = delta*1
            break
        elif F_P < F_Q and F_Q < F_R:
            delta = delta*-1
            break
        else:
            r = random.uniform(a, b)
            delta = random.uniform(10**-9, 10**-12)
            F_P = my_function((r-delta), d, c)
            F_Q = my_function(r, d, c)
            F_R = my_function((r + delta), d, c)

    x_array = [r]

    out.write("#It: "+str(0)+"\t\t"+"r-delta: " + str(round(r - delta, 4))+"\t\t"+"r: " + str(round(r, 4)) + "\t\t"+"r+delta: " + str(round(r +
              delta, 4)) + "\t\t"+"F(r-delta): " + str(round(F_P, 4)) + "\t\t"+"F(r): " + str(round(F_Q, 4)) + "\t\t"+"F(r+delta): " + str(round(F_R, 4)))
    out.write("\n")

    x = r + delta*pow(2, K)  # new value of x
    F_N = my_function(x, d, c)
    x_p = r  # previous value of x
    x_array.append(x)

    out.write("#It \t\t\tx\t\t\tx_p\t\t\tF_N\t\t\tF_Q")
    out.write("\n")

    out.write(str(K)+"\t\t\t" + str(round(x, 4))+"\t\t\t" + str(round(x_p, 4)
                                                                ) + "\t\t\t"+str(round(F_N, 4)) + "\t\t\t"+str(round(F_Q, 4)))
    out.write("\n")

    # main function
    while F_N <= F_Q:
        K = K + 1

        F_Q = F_N
        x_p = x
        x = x + delta*pow(2, K)
        x_array.append(x)
        F_N = my_function(x, d, c)

        out.write(str(K)+"\t\t\t" + str(round(x, 4))+"\t\t\t" + str(round(x_p, 4)
                                                                    ) + "\t\t\t"+str(round(F_N, 4)) + "\t\t\t"+str(round(F_Q, 4)))
        out.write("\n")

    # storing and writing the final boundary values
    b = round(x, 4)
    a = round(x_p-delta*pow(2, K), 4)
    if a > b:
        t = a
        a = b
        b = t
    out.write("\n")
    out.write("Bounding_Phase_iterations are over ")
    out.write("with total " + str(K+3) + " function evaluations.")

    out.write("\n")
    out.write("\n")
    out.write("Now we move to Golden Section Method")
    out.write("\n")
    out.write("\n")
    out.write("Range for the Golden search is: " + str(a)+"," + str(b))
    print("("+str(a)+","+str(b)+")")

    plot_x_versus_iterations(x_array, K, "X", c, "bounding_phase_method")
    plot_range(org_a, org_b, a, b, c, d, "bounding_phase_method")

    return golden_search(a, b, d, c)


# calculate absolute value of x using normalized values of a and b
def calc(w, a, b):
    x = (b-a)*w+a
    return x

# golden search function


def golden_search(a, b, d, c):
    org_a = a
    org_b = b
    a1 = a
    b1 = b
    eps = random.uniform(10**-9, 10**-12)
    aw = 0
    bw = 1
    lw = 1
    w1 = aw+0.618*lw
    w2 = bw-0.618*lw
    # initial function evaluations
    F_w1 = my_function(calc(w1, a1, b1), d, c)
    F_w2 = my_function(calc(w2, a1, b1), d, c)
    K = 1
    t1_array = [round(calc(aw, a1, b1), 4)]
    t2_array = [round(calc(bw, a1, b1), 4)]

    print("************")
    print("Golden Section Search Method")
    out = open(f"Golden_search_iterations_part{c}.out", "w")
    out.write("a: " + str(round(a1, 4))+"\t\t" + "b: "+str(round(b1, 4)) + "\t\t"+"eps: "+str(round(eps, 4)) + "\t\t"+"aw: "+str(round(aw, 4))+"\t\t"+"bw: " + str(round(bw, 4)) +
              "\t\t"+"lw: "+str(round(lw, 4))+"\t\t"+"w1: "+str(round(w1, 4)) + "\t\t"+"w2: "+str(round(w2, 4)) + "\t\t"+"F_w1: "+str(round(F_w1, 4))+"\t\t"+"F_w2: " + str(round(F_w2, 4)))
    out.write("\n")

    out.write("#It \t\taw\t\tbw\t\tlw\t\tw1\t\tw2\t\tF_w1\t\tF_w2")
    out.write("\n")

    # main function
    while abs(lw) > eps:
        if F_w1 > F_w2:
            bw = w1
            lw = bw-aw
            w1 = w2
            w2 = bw-0.618*lw
            F_w1 = F_w2
            F_w2 = my_function(calc(w2, a1, b1), d, c)
        else:
            aw = w2
            lw = bw-aw
            w2 = w1
            w1 = aw+0.618*lw
            F_w2 = F_w1
            F_w1 = my_function(calc(w1, a1, b1), d, c)
        K = K+1
        out.write(str(K)+"\t\t" + str(round(aw, 4))+"\t\t" + str(round(bw, 4)) + "\t\t"+str(round(lw, 4)) + "\t\t" +
                  str(round(w1, 4)) + "\t\t"+str(round(w2, 4)) + "\t\t"+str(round(F_w1, 4)) + "\t\t"+str(round(F_w2, 4)))
        out.write("\n")
        t1 = round(calc(aw, a1, b1), 4)
        t2 = round(calc(bw, a1, b1), 4)
        t1_array.append(t1)
        t2_array.append(t2)

    plot_x_versus_iterations(t1_array, K-1, "t1", c, "golden_search_method")
    plot_x_versus_iterations(t2_array, K-1, "t2", c, "golden_search_method")
    plot_range(org_a, org_b, t1, t2, c, d, "golden_search_method")
    out.write("The final bracketed range is: " +
              str(round(t1))+","+str(round(t2)))
    out.write("The Golden section search completed with " +
              str(K+1)+" iterations.")
    print("("+str(t1)+","+str(t2)+")")
    return (t1, t2)


# main function
def main():
    # asking for inputs
    a = float(input("Enter a = lower limit of x "))
    b = float(input("Enter b = upper limit of x "))
    c = input("Choose the function(1-6): ")
    if c >= "7":
        print("Choose a suitable input b/w 1-6: ")
        main()
    print("Do you want to minimize the function (Y/N)?: ")
    d = input("d = Input your choice: ")
    # calling first function
    bounding_phase(a, b, d, c)
    plt.show()


main()
