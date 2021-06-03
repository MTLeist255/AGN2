# Maybe this will let me read strings/floats and put into seperate arrays
# Maybe this will let me read strings/floats and put into seperate arrays
#splits string according to delimeters
import numpy as np
'''
Let's make a function that can split a string
into list according the given delimeters.
example data: cat;dog:greff,snake/
example delimeters: ,;- /|:
'''
def string_to_splitted_array(data,delimeters):
    #result list
    res = []
    # we will add chars into sub_str until
    # reach a delimeter
    sub_str = ''
    for c in data: #iterate over data char by char
        # if we reached a delimeter, we store the result
        if c in delimeters:
            # avoid empty strings
            if len(sub_str)>0:
                # looks like a valid string.
                res.append(sub_str)
                # reset sub_str to start over
                sub_str = ''
        else:
            # c is not a deilmeter. then it is
            # part of the string.
            sub_str += c
    # there may not be delimeter at end of data.
    # if sub_str is not empty, we should att it to list.
    if len(sub_str)>0:
        res.append(sub_str)
    # result is in res
    return res

# test the function.
delimeters = ',;- /|:'
# read the csv data from console.
csv_string = np.loadtxt('csv_test.txt')
#lets check if working.
splitted_array = string_to_splitted_array(csv_string,delimeters)
print(splitted_array)