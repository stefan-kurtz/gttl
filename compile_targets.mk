%.o: %.cpp
	$(CXX) $(TIME_OPTION) -DCXX_FLAGS=\"'$(CXXFLAGS) $(CFLAGS) $(CPPFLAGS)'\" -c $< -o $@ $(CXXFLAGS) $(CFLAGS) $(CPPFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# In the simplest case, an executable consists of 1 object.
%.x: %.o
	$(LD) $(LDFLAGS) $< -o $@ $(LDLIBS)

