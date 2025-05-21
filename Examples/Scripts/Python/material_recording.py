#!/usr/bin/env python3

import os
import warnings
import argparse

import acts
from acts.examples import (
    UniformVertexGenerator,
    CompositeVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)

import acts.examples.geant4
from acts.examples.odd import getOpenDataDetector

import math

try:
    import acts.examples.geant4.geomodel
except ImportError:
    # geomodel is optional for this script
    pass

u = acts.UnitConstants

_material_recording_executed = False

def runMaterialRecording(
    detectorConstructionFactory,
    outputDir,
    tracksPerEvent=1000000,
    s=None,
):
    global _material_recording_executed
    if _material_recording_executed:
        warnings.warn("Material recording already ran in this process. Expect crashes")
    _material_recording_executed = True

    rnd = RandomNumbers(seed=228)

    # vertexNoWindow=UniformVertexGenerator(
        # mins=acts.Vector4(-7.4, 173.56, 16560, 0),
        # maxs=acts.Vector4(6.3, 346.1, 16560, 0),
    # )
    # weightNoWindow = 346.1 - 173.56
    vertexWindow=UniformVertexGenerator(
        mins=acts.Vector4(-6.3, 78.7, 16630, 0),
        maxs=acts.Vector4(7.5, 108.8, 16630, 0),
    )
    #vertexWindow=UniformVertexGenerator(
    #    mins=acts.Vector4(-20, 30, 16630, 0),
    #    maxs=acts.Vector4(20, 150, 16630, 0),
    #)
    # weightWindow = 173.56 - 75.3
    # weights = [weightNoWindow, weightWindow]
    # vertexGen=CompositeVertexGenerator(
        # [vertexNoWindow, vertexWindow],
        # weights
    # )

    evGen = EventGenerator(
        level=acts.logging.INFO,
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=1),
                vertex=vertexWindow,
                particles=ParametricParticleGenerator(
                    pdg=acts.PdgParticle.eInvalid,
                    charge=0,
                    randomizeCharge=False,
                    mass=0,
                    p=(1 * u.GeV, 1 * u.GeV),
                    thetaMin=0,
                    thetaMax=0,
                    numParticles=tracksPerEvent,
                    etaUniform=False,
                ),
            )
        ],
        outputParticles="particles_initial",
        outputVertices="vertices_initial",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    g4Alg = acts.examples.geant4.Geant4MaterialRecording(
        level=acts.logging.INFO,
        detectorConstructionFactory=detectorConstructionFactory,
        randomNumbers=rnd,
        inputParticles=evGen.config.outputParticles,
        outputMaterialTracks="material-tracks",
    )

    s.addAlgorithm(g4Alg)

    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            maxDistance=500,
            prePostStep=True,
            recalculateTotals=True,
            inputMaterialTracks="material-tracks",
            filePath=os.path.join(outputDir, "geant4_material_tracks_validation.root"),
            level=acts.logging.INFO,
        )
    )

    return s

def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to generate"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Particle tracks per event"
    )
    p.add_argument(
        "-i", "--input", type=str, default="", help="input (GDML/SQL) file (optional)"
    )

    args = p.parse_args()

    detectorConstructionFactory = None
    if args.input.endswith(".gdml"):
        detectorConstructionFactory = (
            acts.examples.geant4.GdmlDetectorConstructionFactory(args.input)
        )
    elif args.input.endswith(".sqlite") or args.input.endswith(".db"):
        geoModelTree = acts.geomodel.readFromDb(args.input)
        detectorConstructionFactory = (
            acts.examples.geant4.geomodel.GeoModelDetectorConstructionFactory(
                geoModelTree
            )
        )

    runMaterialRecording(
        detectorConstructionFactory=detectorConstructionFactory,
        tracksPerEvent=args.tracks,
        outputDir=os.getcwd(),
        s=acts.examples.Sequencer(events=args.events, numThreads=1, trackFpes=False),
    ).run()


if "__main__" == __name__:
    main()
