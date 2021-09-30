/*CITCOM: A finite element convection program written at Caltech 1992 */
/*Aims to include an iterative matrix solver based on Multigrid techniques */
/*To do this requires the use of a mixed method and a conjugate-gradient */
/*approach to determining the */

//
// A particle method is implemented by Shijie Zhong in July 2002.
//
// The first paper that used this particle code from our group
//  is:
// Allen K. McNamara and Shijie Zhong, The influence of thermochemical
//   convection on the fixity of mantle plumes, EPSL, 222, 485-500, 2004
// If you use this code, you may reference this paper and also L. Moresi
// 's original paper on Cartesian citcom.
//
// The particle method works significantly better than the field method
// in my view. Although there is a version of field method implemented
// in this code, I do not recommend to use it.
//                          a note by SZ
#include <mpi.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

main(argc, argv) int argc;
char **argv;

{ /* Functions called by main*/
    void general_stokes_solver();
    void read_instructions();
    void solve_constrained_flow();
    void solve_derived_velocities();
    void process_temp_field();
    void process_heating();

    double dot();

    int k, i, *temp;
    double CPU_time0(), time, initial_time, start_time;

    struct All_variables E;
    /*	parallel_process_initialization(&E,argc,argv); */
    E.parallel.me = 0;
    E.parallel.nproc = 1;
    E.parallel.me_loc[1] = 0;
    E.parallel.me_loc[2] = 0;
    E.parallel.me_loc[3] = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(E.parallel.me));
    MPI_Comm_size(MPI_COMM_WORLD, &(E.parallel.nproc));

    gethostname(E.parallel.machinename, 160);

    E.monitor.solution_cycles = 0;

    if (E.parallel.me == 0)
    {
        start_time = time = CPU_time0();
    }

    read_instructions(&E, argc, argv);

    E.control.keep_going = 1;

    if (E.parallel.me == 0)
    {
        fprintf(stderr, "Input parameters taken from file '%s'\n", argv[1]);
        fprintf(stderr, "Initialization complete after %g seconds\n\n", CPU_time0() - time);
        fflush(E.fp);
        initial_time = CPU_time0() - time;
        fprintf(E.fp, "Initialization overhead = %f\n", initial_time);
        initial_time = CPU_time0();
    }

    general_stokes_solver(&E);
    process_temp_field(&E, E.monitor.solution_cycles);
    process_new_velocity(&E, E.monitor.solution_cycles);

    if (E.control.stokes)
    {
        E.control.keep_going = 0;
        E.monitor.solution_cycles++;
    }

    while (E.control.keep_going && (Emergency_stop == 0))
    {

        E.monitor.solution_cycles++;
        if (E.monitor.solution_cycles > E.control.print_convergence)
            E.control.print_convergence = 1;

        process_heating(&E);

        (E.next_buoyancy_field)(&E);

        process_temp_field(&E, E.monitor.solution_cycles);

        E.monitor.elapsed_time += E.advection.timestep;

        general_stokes_solver(&E);
        process_new_velocity(&E, E.monitor.solution_cycles);

        if (E.control.imposevelo && E.monitor.elapsed_time * E.control.timescale >= E.control.age_total_double)
        {
            E.control.keep_going = 0;
        }

        if (E.control.composition && strcmp(E.control.comp_adv_method, "particle") == 0)
            (E.next_buoyancy_field)(&E);
        for (i = 1; i <= E.mesh.nno; i++)
        {
            E.T_old[i] = E.T[i];
        }

        if (E.parallel.me == 0)
        {
            fprintf(E.fp, "CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e  maxT = %.4e sub_iteration%d markers=%d\n", CPU_time0() - start_time, CPU_time0() - time, E.monitor.solution_cycles, E.monitor.elapsed_time, E.advection.timestep, E.monitor.T_interior, E.advection.last_sub_iterations, E.advection.markers_g);
            time = CPU_time0();
        }
    }

    if (E.parallel.me == 0)
    {
        time = CPU_time0() - initial_time;
        fprintf(E.fp, "Average cpu time taken for velocity step = %f\n", time / ((float)(E.monitor.solution_cycles - 1)));
        fprintf(stderr, "Average cpu time taken for velocity step = %f\n", time / ((float)(E.monitor.solution_cycles - 1)));
    }

    fclose(E.fp);

    parallel_process_termination();

    return;
}
