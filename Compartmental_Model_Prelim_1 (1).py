#!/usr/bin/env python
# coding: utf-8

# # Preliminary Compartmental Model 1 
# still subject to change, troubleshooting, and commentary

# In[1]:


#slider 
#!pip install seir
#!pip install ipywidgets
#!jupyter nbextension enable --py widgetsnbextension
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from IPython.display import display, clear_output, set_matplotlib_formats
set_matplotlib_formats('pdf', 'svg')


# #### Understanding a Compartmental Model 
# 
# $$ \frac{d[D]}{dt} = \zeta [C] - (\alpha [S][D] - \beta [C]),$$
# $$ \frac{d[S]}{dt} = \alpha [S][D],$$
# $$ \frac{d[C]}{dt} = \alpha [S][D] - \beta [C],$$
# $$ \frac{d[R]}{dt} = \beta [C],$$
# $$ \frac{d[T]}{dt} = \gamma [R],$$
# $$ \frac{d[E]}{dt} = 1- \gamma [T],$$
# 
# 
# where $[S]_{max}$ is called the "carrying capapacity"

# More interesting than the compartments of the models are the rules in which they interact:
# $$S + D \to C,$$
# $$C \to R,$$
# $$R \to T + E .$$

# In[2]:


##slider goes here
def graph(alpha=1, beta=1, gamma=1, zeta=1):
# SIR Model

    def dSIcRd(current_state, time):
    
        donor_cell = current_state[0]
        susceptible_cell = current_state[1]
        conjugated_cell = current_state[2]
        recipient_cell = current_state[3]
        transformed_cell = current_state[4]
        dead_cell = current_state[5]

        susceptible_cell_dot = -alpha*susceptible_cell*donor_cell
        conjugated_cell_dot = alpha*susceptible_cell*donor_cell - beta*conjugated_cell
        recipient_cell_dot = beta*conjugated_cell
        transformed_cell_dot = gamma*recipient_cell
        dead_cell_dot = (1 - gamma)*recipient_cell
        donor_cell_dot = zeta*conjugated_cell - (alpha*susceptible_cell*donor_cell - beta*conjugated_cell)
        #look into zeta 


        return [donor_cell_dot, susceptible_cell_dot, conjugated_cell_dot, 
            recipient_cell_dot, transformed_cell_dot, dead_cell_dot]

    

#Consider Initial Conditions
    initial_conditions = [100,1000,0,0,0,0]
    time = np.linspace(0,20,1000)

    result = odeint(dSIcRd, initial_conditions, time)

    donor = result[:,0]
    susceptible = result[:,1]
    conjugated = result[:,2]
    recipient = result[:,3]
    transformed = result[:,4]
    dead = result[:,5]

    plt.plot(time, donor, label='Donor Cells')
    plt.plot(time, susceptible, label='Susceptible Cells')
    plt.plot(time, conjugated, label='Cells Conjugated')
    plt.plot(time, recipient, label='Cells Recieving Donor Gene')
    plt.plot(time, transformed, label='Cells Transformed')
    plt.plot(time, dead, label='Cell Death')


    plt.grid()
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Population (thousands)')
    plt.title('Solution of the dSIcRd Model')
interact(graph,alpha=(0.000,1,0.001), beta=(0.000,1,0.001), gamma=(0.000,1,0.001), zeta=(0.000,1,0.001))


# In[ ]:




