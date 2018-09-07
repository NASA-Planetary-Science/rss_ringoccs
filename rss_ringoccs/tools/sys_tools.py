import subprocess
import os
import time

def shell_execute(script):
    """
        Function:
            shell_execute
        Purpose:
            Execute a shell script from within Python.
        Variables:
            script:     A list containing the path to the shell
                        script and variables. For example:
                            script = ['path/to/script','v1',...,'vn']
                        Elements must be strings.
        Outputs:
            Process:    An instance of the Popen class from the
                        subprocess module. This contains attributes
                        such as the arguments passed to it, and other
                        system details.
        Dependencies:
            [1] subprocess
        Notes:
            This routine has only been tested using scripts written
            in Bash, and on standard Unix commands.
        References:
            [1] https://docs.python.org/3/library/subprocess.html
            [2] https://stackoverflow.com/questions/
                3777301/how-to-call-a-shell-script-from-python-code
        Examples:
            Suppose we have the following shell script test.sh:
                #!/bin/bash
                printf "Hello, World! My name is %s!" "$1"
            Run this shell script inside of Python:
                In [1]: import diffcorr as dc
                In [2]: dc.shell_execute(['./test.sh','Bob'])
                        Hello World! My name is Bob!
            We can also execute simple Unix commands.
                In [1]: import diffcorr as dc
                In [2]: a = dc.shell_execute(["echo","Bob"])
                        Bob
                In [3]: a.args
                Out[3]: [' echo Bob']
        History:
            Created: RJM - 2018/05/16 5:49 P.M.
    """
    string = ""
    for x in script:
        if (not isinstance(x, str)):
            raise TypeError("Input must be a list of strings")
        else:
            string = "%s %s" % (string,x)
    Process=subprocess.Popen([string],shell=True)
    return Process

def date_string():
    """
        Function:
            date_string
        Purpose:
            Create string of the form "yyyy_mm_dd_HH_MM_SS_"
        Variables:  There are no variables to this function.
        Outputs:
            date:   Current date "year/month/day/hour/minute/second"
        Dependencies:
            [1] time
        Notes:
            The end of the string has an underscore "_"
        Examples:
            Get the current date.
                In [1]: import diffcorr as dc
                In [2]: dc.date_string()
                Out[2]: '2018_05_27_11_28_22_'
        History:
            Created: RJM - 2018/05/27 11:27 A.M.
    """
    strings = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = strings.split(',')
    date=""
    for x in t:
        date = "%s%s_" % (date,x)
    return date

def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)
