/*
To choose a color that is somewhat darker or lighter than another
color, the function |choose_color| queries the RGB values of a pixel
and multiplies them with a factor. If all goes well, the function
returns |True|. If the chosen color ends up being the same as the
original, the color gray75 is returned instead.
*/

#include <stdio.h>
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include "Converters.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

Boolean choose_color(self, factor, base, result)
Widget self;
double factor;
Pixel base;
Pixel *result;
{
    Colormap colormap;
    XColor color, dummy;

    if (XtIsRealized(self))
	colormap = self->core.colormap;
    else
	colormap = DefaultColormapOfScreen(XtScreen(self));
    color.pixel = base;

    XQueryColor(XtDisplay(self), colormap, &color);
    color.red = min(65535, factor * color.red);
    color.green = min(65535, factor * color.green);
    color.blue = min(65535, factor * color.blue);
    if (! XAllocColor(XtDisplay(self), colormap, &color))
	return False;
    if (base == color.pixel
	&& ! XAllocNamedColor(XtDisplay(self), colormap, "gray75",
			      &color, &dummy))
	return False;

    *result = color.pixel;
    return True;

}

