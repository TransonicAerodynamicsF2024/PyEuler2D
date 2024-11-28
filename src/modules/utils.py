import platform
import csv
import matplotlib.pyplot as plt

def console_printer(iteration, Cl, Cd, Cm):
    if not hasattr(console_printer, "header_printed"):
        console_printer.header_printed = False

    if iteration % 50 == 0 or not console_printer.header_printed:
        print("\n{:<10} {:>15} {:>15} {:>15}".format("Iteration", "Cl", "Cd", "Cm"))
        print("-" * 60)
        console_printer.header_printed = True

    print("{:<10d} {:>15.6f} {:>15.6f} {:>15.6f}".format(iteration, Cl, Cd, Cm))

def print_green(text):
    """
    Function to print a string in Green color
    Args:
        text (str): String to colorize green
    """
    attr = ['42', '1']
    if platform.system() == "Windows":
        print(text)
    else:
        print('\x1b[%sm%s\x1b[0m' % (';'.join(attr), text))


def text_green(text):
    """
    Function to print a string in Green color
    Args:
        text (str): String to colorize green
    """
    attr = ['42', '1']
    if platform.system() == "Windows":
        return text
    else:
        return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), text)


def print_red(text):
    """
    Function to print a text in Red for error notifications
    Args:
        text (str): String to colorize in Red
    """
    attr = ['41', '1']

    if platform.system() == "Windows":
        print(text)
    else:
        print('\x1b[%sm%s\x1b[0m' % (';'.join(attr), text))


def text_red(text):
    """
    Function to print a text in Red for error notifications
    Args:
        text (str): String to colorize in Red
    """
    attr = ['41', '1']

    if platform.system() == "Windows":
        return text
    else:
        return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), text)


def printHead():
    head="""
    .         .                   .                                                      
       .                                                                                     
               .   .                   .                                                 
    PyEuler 2D 
    ------                      .       
    PROGRAMMED FOR TRANSONIC AERODYNAMICS FALL 2024                           .
           .                 .                               .                            
      
                         .                .                     . .                      
                    .    .         .                                                     
                                                                                        .
                  .                      .      .                                       .
                   .                                  .                                  
    """

    print(head)
    print("-"*100)

