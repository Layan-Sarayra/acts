#!/usr/bin/env python3
import os, argparse, pathlib, acts, acts.examples
from pathlib import Path
from typing import Optional, Union

parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector")
parser.add_argument("--events", "-n", help="Number of events", type=int, default=10)
args = vars(parser.parse_args())

outputDir = pathlib.Path("/eos/user/l/lalsaray/KDE_output")

s = acts.examples.Sequencer(events=args["events"], numThreads=1, outputDir=str(outputDir))

@acts.examples.NamedTypeArgs(logLevel=acts.logging.Level)
def KDE_printer(
    s,
    logLevel: Optional[acts.logging.Level] = None,
    # outputDirRoot: Optional[Union[Path, str]] = None,   
) -> None:
    from acts.examples import KDEAlgorithm
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    KDE = KDEAlgorithm(level = customLogLevel())
    s.addAlgorithm(KDE)

    return s

if __name__ == "__main__":
    args = vars(parser.parse_args())
    os.environ["ACTS_SEQUENCER_DISABLE_FPEMON"] = "1"

    s = acts.examples.Sequencer(events=args["events"], numThreads=1, outputDir=str(outputDir))

    s = KDE_printer(s)  # Add KDE printing to the sequencer
    
    s.run()    # Execute the sequencer
