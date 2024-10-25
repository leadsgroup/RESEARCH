class Example_class():
    """
    Description of the class

    Attributes
    ----------
    Describe non-method attributes here
    x : data_type
        Description of the attribute
    a : data_type (float, int, str, etc)
        Description of the attribute a

    Methods
    -------
    Only list public methods here that can and will be called by the user.
    Do not list secondary (aka private) methods that are only called by other functions.
    
    method1
        Description of the return.

    Notes
    -----
    Additional information that is not strictly necessary for the user, but
    may be helpful. 
    
    **Definitions**
    'term'
        term definition
    
    References
    ----------

    .. [1] 
    
    """
    
    def  method1(self):
        print("Hello World")
        self.private_method2()
        
    def private_method2(self):
        # This is private
        pass    