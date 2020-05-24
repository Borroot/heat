all:
	@$(MAKE) --silent -C seq
	@$(MAKE) --silent -C openmp
	@$(MAKE) --silent -C mpi

clean:
	@$(MAKE) clean --silent -C seq
	@$(MAKE) clean --silent -C openmp
	@$(MAKE) clean --silent -C mpi

.PHONY: all clean
