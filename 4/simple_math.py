def add(a,b):
    """
    Takes two numbers and returns the sum. 
    
    Parameters
    ----------
    a : float
        First number to add. 
    b : float
        Second number to add.

    Returns
    -------
    c : float
        The c of a and b.
    """

    c = a + b
    return c

def square_add(a,b):
    """
    Takes two numbers and add their squares. 
    
    Parameters
    ----------
    a : float
        First number to square and add. 
    b : float
        Second number to square and add.

    Returns
    -------
    c : float
        The resulting value.
    """
    c = a**2 + b**2
    return c

def quotrem(a,b):
    """
    Takes two numbers and returns the quotient and remainder from the division.

    Parameters
    ----------
    a : float
        Numerator in the division.
    b : float
        Denominator in the division.
    
    Returns
    -------
    quot : float
        Quotient from the division.
    rem : float
        Remainder from the division.
    """
    quot = a % b
    rem = a // b
    return quot, rem

