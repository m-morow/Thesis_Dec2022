## Madeleine Tan 
#### mmtan@umich.edu
#### Undergraduate Thesis Dec. 2022
--------------------------------------- 
##### Code:
        slant_stack.py
                calculates simple linear moveout correction based on distance to event
                returns: figure with 3 subplots [record section, stack, scatter plot of distance of event]        
        plot_rf.py
                plots record section of receiver functions
                returns: figure with 3 subplots [record section, stack, scatter plot of back-azimuth of event]
        hk.py
                calculates H and kappa values of receiver function stream
                returns: contourf figure of H (x axis) and kappa (y axis) and optimum H, k values
        stacking_by_slowness:
                stack traces at each station by slowness
                returns: figure with 5 stacks; all slowness values in s/km, [0.04-0.05), [0.05-0.06), [0.06-0.07), [0.07-0.08]
        azeq_map:
                plot events on azimuthal equidistant map projection
                returns: figure with event locations
                

##### Reference code:
        SeisPy by Mijian Xu: https://github.com/xumi1993/seispy
        
##### Figures:
        Record section by station ./rf_figures
        Hk ./figures
        diagrams ./diagrams
        
        
