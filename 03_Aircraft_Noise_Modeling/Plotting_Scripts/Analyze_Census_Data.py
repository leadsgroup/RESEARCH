

# Python imports
import matplotlib.pyplot as plt  
import sys 
import numpy as np 
import json 
 
def main(): 
    
    # Open and read the JSON file
    with open('data.json', 'r') as file:
        data = json.load(file)
    
    # Print the data
    print(data)
    
    return


if __name__ == '__main__': 
    main()    
    plt.show()
