# windows setup
CXX = g++
CXXFLAGS = -g
EXEC = main
OBJDIR=dist1
OBJS = $(addprefix $(OBJDIR)/, \
	main.o)

${EXEC}: ${OBJS}
	${CXX} ${CXXFLAGS} -o $(OBJDIR)/${EXEC} ${OBJS}


$(OBJDIR)/%.o: %.cpp
	if not exist $(OBJDIR) mkdir $(OBJDIR)
	${CXX} ${CXXFLAGS} -o $@ -c $<

.PHONY: clean

run:
	$(OBJDIR)/${EXEC}

clean:
	del /S /Q $(OBJDIR)
	rmdir $(OBJDIR)