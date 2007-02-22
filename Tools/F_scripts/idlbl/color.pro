function color, name, GET_LIST = color_list

; this function returns the color index corresponding to the input name
;
; We reserve the first 12 colors of the colormap.  The names of
; these colors are stored in this function.  This allows these
; routines to use the same color for specific tasks (eg. background
; color),  regardless of the colormap used in the plotting.
;
; optional argument: The name of available colors can be returned
; through the GET_LIST argument.
;

; define the standard colors
std_colors = ['black', $
              'ltgreen', $
              'dkgreen', $
              'dkblue', $
              'ltblue', $
              'yellow', $
              'red', $
              'dkgray', $
              'mdgray', $
              'ltgray', $
              'verygray', $
              'white']

color_index = (where(std_colors EQ name))[0] > 0

color_list = std_colors

return, color_index
end

