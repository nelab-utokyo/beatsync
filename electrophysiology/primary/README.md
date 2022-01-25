
# data

mean MUA of 5-30 ms after stimuli.
・evoked_res2 for periodic and rhythmic 
・evoked_resm2 for music condition

mean MUA of 5-30 ms subtracted by 0-4ms after stimuli.
・evoked_res3 for periodic and rhythmic 
・evoked_resm3 for music condition

responsep.mat
・mean MUA response for on beat and off beat in periodic condition
・rmv is the TC electrode without click response

responsef.mat
・mean MUA response for on beat and off beat in rhythmic condition

responsem.mat
・mean MUA response for on beat and off beat in music condition
・rmv is the TC electrode with longer spike latency (> 20ms)


# code

## regionmap.m

main function for determining A1/AAF/belt region.
Variables, "corebelt" & "a1aaf", are added to CF_latency.mat
<br>
(corebelt == 1 for core, == 0 for belt & non auditory)
<br>
(a1aaf == 1 for a1, == -1 for aaf, 0 for belt & non auditory)

## set_threshold.m 

main function for processing individual data.
<br>
evoked_res2 is used for thresholding MUA with click response.
<br>
evoke_res3, evoke_resm3 -> responsep, responsef, responsem
<br>
Also generates the nessecary variable for fitting, experimental_data.mat

## alldata_plot.m

main function for plotting the whole individual data.
<br>
responsep, responsef, responsem -> plot (Fig. 3c-e, Fig. 4c)

## contrast_map.m

analyze the beat contrast depending on the region.
Only Core/belt is categorized.
<br>
responsep, responsef, responsem -> plot (Extended Data Fig. 3b)

## individual_plot.m

Plot the individual data.
<br>
responsep, responsef, responsem -> plot (Fig. 3b, Fig. 4b, Extended Data Fig. 3a)

## random_rawplot.m

Plot the MUA data of random click.
(Fig. 4h)


