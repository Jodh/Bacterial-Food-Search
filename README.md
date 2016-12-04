# Bacterial-Food-Search
A distributed gradient decent algorithm to mimic bacterial food search, written as a project for COT5405. 
My first attempts at actuall programming, hope all goes well.
The original code can be found here: https://github.com/sss1/bact-sim , but it is in matlab and incomplete and the functions are incorrectly implemented. In short the code is out-dated. It looks like they abandoned it in favour of locally developing the code.

If you are looking to implement it yourself please use the original research paper by the same people, which can be found here: http://www.snl.salk.edu/~navlakha/pubs/recomb2016.pdf

# Running and understanding the output
1. Download DGD.py
2. Install numpy, matplotlib and scipy libraries
  1. For Ubuntu use command ```sudo apt-get install python-name_of_library```
  2. Replace *name_of_library* by numpy, matplotlib and scipy respectively to install them
3. Open terminal, change directory to where ever DGD.py is stored
4. Enter this in the terminal ```python DGD.py``` 
5. Three .png image files will be generated
  1. before.png : The state of the swarm in the beginning
  2. after.png  : State of the swarm after the simulation has run
  3. timeline.png : Gives the number of bacteria that have reached the food source per iteration
