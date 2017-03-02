#Compiler to use: gcc version 2.96 20000731 (Red Hat Linux 7.1 2.96-85)
CCC=/usr/bin/g++
CC=/usr/bin/gcc

#Application name
NAME = dragon

#directories
XmGraph = XmGraph
Xfwf = Xfwf
Main_Dragon = Main_Dragon
apolist_dir = Main_Dragon/apolist
util_dir    = Main_Dragon/util
vcg_dir = vcg

#Source Files
SRCFILES = $(Main_Dragon)/UIClasses.cxx \
	$(Main_Dragon)/cfgwindow.cxx  \
	$(Main_Dragon)/cfgwindow_stubs.cxx \
	$(Main_Dragon)/graphwindow.cxx  \
	$(Main_Dragon)/dragon_stubs.cxx  \
	$(Main_Dragon)/graphwindow_stubs.cxx  \
	$(Main_Dragon)/dragon.cxx	\
	$(Main_Dragon)/FGnode.cxx	\
	$(Main_Dragon)/ddgui.cxx        \
	$(Main_Dragon)/regions.cxx      \
	$(apolist_dir)/apogen.cxx      \
	$(apolist_dir)/apogui.cxx      \
	$(apolist_dir)/Apolist.cxx      \
	$(util_dir)/readfile.cxx      \
	$(vcg_dir)/GraphBasic.cxx	\
	$(vcg_dir)/GraphEdge.cxx	\
	$(vcg_dir)/GraphNode.cxx	\
	$(vcg_dir)/GraphVCG.cxx		\
	$(vcg_dir)/PrintGraph.cxx		


XPACKSRC = $(XmGraph)/Graph.c $(XmGraph)/Arc.c $(Xfwf)/Board.c $(Xfwf)/DrawImageString.c $(Xfwf)/Frame.c $(Xfwf)/Slider2.c $(Xfwf)/TextWidth.c $(Xfwf)/scroll.c $(Xfwf)/strnchr.c $(Xfwf)/Common.c $(Xfwf)/DrawString.c $(Xfwf)/Label.c $(Xfwf)/Tablist2Tabs.c $(Xfwf)/choosecol.c $(Xfwf)/strarray.c

OBJECTS =  $(XPACKSRC:.c=.o) $(SRCFILES:.cxx=.o)
DEPFILES = $(OBJECTS:.o=.d)

LIBRARY = $(vcg_dir)/libPrintGraph.so
LIBRARY =


OPT = -D_XFWFINTERFACE -DNO_TIFF -DNO_JPEG -DNO_RLE -DNO_GIF -DNO_PBM -DNO_XPM -D__USE_FIXED_PROTOTYPES__ -DPP_GUI -DCT_PARAMETER -DTRY_SPRING -DLINUX -DXM_1_1_BC

DEBUG = -g

#GUI include directories and libraries to link
GUIINCL = -I. -I/usr/X11R6/LessTif/Motif1.2/include -I/usr/X11R6/include/X11 -IXfwf  -I/usr/include/mysql -I /usr/local/mysql/include/mysql/ -I/usr/X11R6/include
GUILNK =  -L/usr/X11R6/LessTif/Motif1.2/lib -lXm -L/usr/X11R6/lib -lXt -lXbae -lX11 -lXext -lXmu -lm  -L /usr/lib/mysql -L /usr/local/mysql/lib/mysql -lmysqlclient -Wl,-R/usr/local/mysql/lib/mysql

#Linking Part
$(NAME): $(XPACKSRC) $(SRCFILES) $(OBJECTS) 
	$(CCC) $(GUILNK) $(DEBUG) -o $(NAME) $(OBJECTS) $(LIBRARY)


#Compiling Part
%.o: %.cxx
	$(CCC) -c $(OPT) $(GUIINCL) $(DEBUG) $< -o $@

%.o: %.c
	$(CC) -c $(OPT) $(GUIINCL)  $(DEBUG) $< -o $@
 
#Clean Part
clean:
	rm -f .*~ core $(NAME) $(OBJECTS) $(DEPFILES)
sclean:
	rm -f core $(Main_Dragon)/*.o


#Dragon self analyzer
self: $(XPACKSRC) $(SRCFILES) $(OBJECTS) 
	opencc -ipa -dragon -O2 $(GUILNK) $(DEBUG) -o dragon.self $(OBJECTS) $(LIBRARY)




