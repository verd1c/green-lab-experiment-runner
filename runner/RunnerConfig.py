import os
import signal
from EventManager.Models.RunnerEvents import RunnerEvents
from EventManager.EventSubscriptionController import EventSubscriptionController
from ConfigValidator.Config.Models.RunTableModel import RunTableModel
from ConfigValidator.Config.Models.FactorModel import FactorModel
from ConfigValidator.Config.Models.RunnerContext import RunnerContext
from ConfigValidator.Config.Models.OperationType import OperationType
from ProgressManager.Output.OutputProcedure import OutputProcedure as output

from typing import Dict, List, Any, Optional
from pathlib import Path
from os.path import dirname, realpath

import pandas as pd
import time
import subprocess
import shlex


class RunnerConfig:
    ROOT_DIR = Path(dirname(realpath(__file__)))

    # ================================ USER SPECIFIC CONFIG ================================
    """The name of the experiment."""
    name:                       str             = "new_runner_experiment"

    """The path in which Experiment Runner will create a folder with the name `self.name`, in order to store the
    results from this experiment. (Path does not need to exist - it will be created if necessary.)
    Output path defaults to the config file's path, inside the folder 'experiments'"""
    results_output_path:        Path             = ROOT_DIR / 'experiments'

    """Experiment operation type. Unless you manually want to initiate each run, use `OperationType.AUTO`."""
    operation_type:             OperationType   = OperationType.AUTO

    """The time Experiment Runner will wait after a run completes.
    This can be essential to accommodate for cooldown periods on some systems."""
    time_between_runs_in_ms:    int             = 1000

    # some static config values
    DIR_SYSTEMS: Path = Path("/home/leafyit/tmp/jason/repo/systems")
    SYSTEM_COARSE_GRAINED_PATH = "/home/leafyit/tmp/jason/repo/systems/new"
    SYSTEM_ALL_ATOM_PATH = "/home/leafyit/tmp/jason/repo/systems/1ka_aa"
    FILE_TMP_RUN_DIR: Path = ROOT_DIR / 'tmp_run_files'

    # Dynamic configurations can be one-time satisfied here before the program takes the config as-is
    # e.g. Setting some variable based on some criteria
    def __init__(self):
        """Executes immediately after program start, on config load"""

        EventSubscriptionController.subscribe_to_multiple_events([
            (RunnerEvents.BEFORE_EXPERIMENT, self.before_experiment),
            (RunnerEvents.BEFORE_RUN       , self.before_run       ),
            (RunnerEvents.START_RUN        , self.start_run        ),
            (RunnerEvents.START_MEASUREMENT, self.start_measurement),
            (RunnerEvents.INTERACT         , self.interact         ),
            (RunnerEvents.STOP_MEASUREMENT , self.stop_measurement ),
            (RunnerEvents.STOP_RUN         , self.stop_run         ),
            (RunnerEvents.POPULATE_RUN_DATA, self.populate_run_data),
            (RunnerEvents.AFTER_EXPERIMENT , self.after_experiment )
        ])
        self.run_table_model = None  # Initialized later
        output.console_log("Custom config loaded")

    def create_run_table_model(self) -> RunTableModel:
        """Create and return the run_table model here. A run_table is a List (rows) of tuples (columns),
        representing each run performed"""

        #cpu_limit_factor = FactorModel("force_field", [ "allatom", "coarsegrained" ])
        #alg_factor  = FactorModel("algorithm" , [ "leapfrog", "steepdecent" ])
        #em_steps_factor  = FactorModel("em_steps" , [ 500, 2500, 5000 ])
        cpu_limit_factor = FactorModel("force_field", [ "charmm", "allatom" ])
        alg_factor  = FactorModel("algorithm" , [ "md", "sd" ])
        md_steps_factor  = FactorModel("md_steps" , [ 50000, 250000, 500000 ])
        self.run_table_model = RunTableModel(
            factors = [cpu_limit_factor, alg_factor, md_steps_factor],
            repetitions=10,
            # exclude_variations = [
            #     {cpu_limit_factor: [70], pin_core_factor: [False]} # all runs having the combination <'70', 'False'> will be excluded
            # ],
            data_columns=[ 'exec_time', 'avg_cpu', 'avg_mem', 'energy_usage', 'accuracy' ]
        )
        return self.run_table_model

    def before_experiment(self) -> None:
        """Perform any activity required before starting the experiment here
        Invoked only once during the lifetime of the program."""
        
        self.FILE_TMP_RUN_DIR.mkdir(parents=True, exist_ok=True)
        

    def before_run(self) -> None:
        """Perform any activity required before starting a run.
        No context is available here as the run is not yet active (BEFORE RUN)"""

        # we want to wipe the FILE_TMP_RUN_DIR before the experiment
        for root, dirs, files in os.walk(self.FILE_TMP_RUN_DIR, topdown=False):
            for file in files:
                os.remove(os.path.join(root, file))

            for dirz in dirs:
                os.rmdir(os.path.join(root, dirz))

        return

    def _get_str_file_for_cfg(self, cfg, stage):
        return f"{stage}_{cfg['force_field']}_{cfg['algorithm']}_{cfg['em_steps']}"

    def _get_alg_path_for_cfg(self, cfg):
        return f"{cfg['algorithm']}_{cfg['md_steps']}"

    def _do_run(self, context: RunnerContext, cfg) -> None:
        """
        http://www.mdtutorials.com/gmx/lysozyme/05_EM.html

        We need to first run:
        gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

        - md.mdp is in algorithms/
        - npt.gro, npt.cpt, topol.top are in the protein dir
        - output md.tpr will just be in tmp folder

        then gmx mdrun -deffnm md_0_1

        """

        path_protein_system_gro = Path(cfg["tgt_sys_path"]) / "1AKI_solv_ions.gro"
        path_protein_system_top = Path(cfg["tgt_sys_path"]) / "topol.top"
        path_protein_system_cpt = Path(cfg["tgt_sys_path"]) / "npt.cpt"
        
        mdp_file_basename = self._get_alg_path_for_cfg(cfg)
        path_alg_mdp = self.DIR_SYSTEMS / "algorithms" / (mdp_file_basename + ".mdp")

        grompp_params = [
            "gmx",
            "grompp",
            "-f",
            str(path_alg_mdp),
            "-c",
            str(path_protein_system_gro),
            "-t",
            str(path_protein_system_cpt),
            "-p",
            str(path_protein_system_top),
            "-o",
            str(self.FILE_TMP_RUN_DIR / ("md_0_1.tpr"))
        ]

        print(f"Running grompp with {grompp_params} rootdir is {str(self.ROOT_DIR)}")
        print(' '.join(grompp_params))

        self.target = subprocess.Popen(
            grompp_params,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.FILE_TMP_RUN_DIR
        )
        out, err = self.target.communicate()
        print(out.decode())
        print(err.decode())


        md_run_params = [
            "gmx",
            "mdrun",
            "-v",
            "-deffnm",
            "md_0_1"
        ]
        print(f"running for params {md_run_params}")
        print(' '.join(md_run_params))
        self.target = subprocess.Popen(
            md_run_params,
                cwd=self.FILE_TMP_RUN_DIR
        )
        # out, err = self.target.communicate()
        # print(out.decode())
        # print(err.decode())

        return

    def start_run(self, context: RunnerContext) -> None:
        """Perform any activity required for starting the run here.
        For example, starting the target system to measure.
        Activities after starting the run should also be performed here."""


        force_field  = context.run_variation['force_field']
        algorithm = context.run_variation['algorithm']
        md_steps = context.run_variation['md_steps']

        target_system_path = self.SYSTEM_ALL_ATOM_PATH if force_field == "allatom" else self.SYSTEM_COARSE_GRAINED_PATH
        print(f'running for {force_field}, {algorithm}, {md_steps}, {target_system_path}')
        run_cfg = {
                "force_field": force_field,
                "algorithm": algorithm,
                "md_steps": md_steps,
                "tgt_sys_path": target_system_path
        }

        self.start_time = time.time()
        em_res = self._do_run(context, run_cfg)

        # cpu_limit = context.run_variation['cpu_limit']
        # pin_core  = context.run_variation['pin_core']

        # start the target
        # self.target = subprocess.Popen(['./primer'],
        #     stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.ROOT_DIR,
        # )

        # Configure the environment based on the current variation
        # if pin_core:
        #     subprocess.check_call(shlex.split(f'taskset -cp 0  {self.target.pid}'))
        # subprocess.check_call(shlex.split(f'cpulimit -b -p {self.target.pid} --limit {cpu_limit}'))
        pass
        

    def start_measurement(self, context: RunnerContext) -> None:
        """Perform any activity required for starting measurements."""

        # man 1 ps
        # %cpu:
        #   cpu utilization of the process in "##.#" format.  Currently, it is the CPU time used
        #   divided by the time the process has been running (cputime/realtime ratio), expressed
        #   as a percentage.  It will not add up to 100% unless you are lucky.  (alias pcpu).
        profiler_cmd = f'ps -p {self.target.pid} --noheader -o %cpu,rss'
        # profiler_cmd = f'ps -p {self.target.pid}'
        wrapper_script = f'''
        while true; do {profiler_cmd}; sleep 1; done
        '''

        time.sleep(1) # allow the process to run a little before measuring
        self.profiler = subprocess.Popen(
            ['sh', '-c', wrapper_script],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        self.pj_profiler = subprocess.Popen(
            shlex.split(
                f'powerjoular -l -p {self.target.pid} -f {context.run_dir / "powerjoular.csv"}'
            )
        )

        pass

    def interact(self, context: RunnerContext) -> None:
        """Perform any interaction with the running target system here, or block here until the target finishes."""

        # No interaction. We just run it for XX seconds.
        # Another example would be to wait for the target to finish, e.g. via `self.target.wait()`
        output.console_log("Running program until it exits seconds")
        self.target.wait()
        # time.sleep(20)

    def stop_measurement(self, context: RunnerContext) -> None:
        """Perform any activity here required for stopping measurements."""

        # out, err = self.profiler.communicate()
        # print(out.decode())
        # print(err.decode())

        self.profiler.kill()
        self.profiler.wait()

        os.kill(self.pj_profiler.pid, signal.SIGINT) # graceful shutdown of powerjoular
        self.pj_profiler.wait()
        pass

    def stop_run(self, context: RunnerContext) -> None:
        """Perform any activity here required for stopping the run.
        Activities after stopping the run should also be performed here."""
        
        self.end_time = time.time()
        self.target.kill()
        self.target.wait()
        pass
    
    def populate_run_data(self, context: RunnerContext) -> Optional[Dict[str, Any]]:
        """Parse and process any measurement data here.
        You can also store the raw measurement data under `context.run_dir`
        Returns a dictionary with keys `self.run_table_model.data_columns` and their values populated"""

        ps_lines = self.profiler.stdout.readlines()
        print(f'time taken: {self.end_time - self.start_time}')

        df = pd.DataFrame(columns=['cpu_usage', 'mem_usage'])
        for i, l in enumerate(ps_lines):
            line = l.decode('ascii').split(' ')
            cpu_usage=float(line[0].strip())
            mem_usage=float(line[1].strip())
            df.loc[i] = [ cpu_usage, mem_usage ]
        
        df.to_csv(context.run_dir / 'raw_data.csv', index=False)

        # get accuracy
        # gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
        p = subprocess.Popen(
            shlex.split(
                'gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center'
            ),
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.FILE_TMP_RUN_DIR
        )
        out, err = p.communicate(input=b'1\n0\n')
        # print(out, err)

        p = subprocess.Popen(
            shlex.split(
                'gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns'
            ),
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.FILE_TMP_RUN_DIR
        )
        out, err = p.communicate(input=b'1\n0\n')
        # print(out, err)

        p = subprocess.Popen(
            shlex.split(
                'gmx analyze -f rmsd.xvg'
            ),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.FILE_TMP_RUN_DIR
        )
        out, err = p.communicate()
        out = out.decode('ascii')

        run_data = {
            'exec_time': round(self.end_time - self.start_time, 3),
            'avg_cpu': round(df['cpu_usage'].mean(), 3),
            'avg_mem': round(df['mem_usage'].mean(), 3)
        }

        for line in out.split("\n"):
            if "SS1" in line:
                sp = line.split(" ")
                for c in sp:
                    if c != "SS1" and c != "":
                        run_data["accuracy"] = c
                        break

        df = pd.read_csv(context.run_dir / f"powerjoular.csv-{self.target.pid}.csv")
        run_data["energy_usage"] = round(df['CPU Power'].sum(), 3)

        force_field  = context.run_variation['force_field']
        algorithm = context.run_variation['algorithm']
        md_steps = context.run_variation['md_steps']
        with open(f"/home/leafyit/tmp/jason/backup.txt", "a+") as f:
            f.write(f"{context.run_nr}_{force_field}_{algorithm}_{md_steps} => {str(run_data)}\n")
            f.close()

        return run_data

    def after_experiment(self) -> None:
        """Perform any activity required after stopping the experiment here
        Invoked only once during the lifetime of the program."""
        pass

    # ================================ DO NOT ALTER BELOW THIS LINE ================================
    experiment_path:            Path             = None
