class ClassName(Parent):
    """
    Brief description of the class

    Attributes
    ----------
    simple_attribute : data_type
        Description of simple attribute
    
    complex_attribute : Data
        Collection of related measurements
        - sub_attribute_1 : data_type
            Description of first sub-attribute
        - sub_attribute_2 : data_type
            Description of second sub-attribute
    
    container_attribute : Container
        Description of container purpose

    Methods
    -------
    method_name(param1, param2)
        Brief description of what the method does

    Notes
    -----
    Additional important information about the class.
    Implementation details, usage guidelines, etc.
    
    **Definitions**
    'Technical Term'
        Definition of the technical term
    'Another Term'
        Definition of another technical term
    
    References
    ----------
    .. [1] Author Name, "Paper Title", Journal, Year
    .. [2] Documentation Link, etc.
    """
    
    def  method1(self):
        print("Hello World")
        self.private_method2()
        
    def private_method2(self):
        # This is private
        pass    