#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
)

u = acts.UnitConstants

def estimateLookup():
    # Initialize the geometry
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[4, 10],
        positions=[30, 60, 90],
        stereos=[0, 0, 0],
        binValue=2,
        surfaceType=0
    )

    # Set up the track lookup grid accumulator
    gridAccumulatorConfig = acts.examples.TrackLookupGridAccumulator.Config(
        bins=[1, 100],
        xBounds=[-1,1],
        yBounds=[-10,10])
    gridAccumulator = acts.examples.TrackLookupGridAccumulator(gridAccumulatorConfig)

    # Set up the track lookup grid writer
    jsonWriterConfig = acts.examples.JsonTrackLookupGridWriter.Config(
        path="/home/romanurmanov/tools/acts/acts_telescope_geo/ActsTelescopeGeometryDevelopment_build/test"
    )
    jsonWriter = acts.examples.JsonTrackLookupGridWriter(jsonWriterConfig)

    # Set up the track estimation algorithm
    trackEstConfig = acts.examples.TrackLookupEstimationAlgorithm.Config(
        refLayers=[list(trackingGeometry.geoIdSurfaceMap().values())[0]],
        inputHits="simhits",
        inputParticles="particles_input",
        outputIPTrackParameters = "IPTrackParameters",
        outputRefLayerTrackParameters = "RefTrackParameters",
        trackLookupGridAccumulator = gridAccumulator,
        trackLookupGridWriters = [jsonWriter]
    )
    trackEstAlg = acts.examples.TrackLookupEstimationAlgorithm(trackEstConfig, acts.logging.INFO)

    # Set up the magnetic field
    field = acts.ConstantBField(acts.Vector3(50 * u.T, 0, 0))

    outputDir = Path.cwd() / "telescope_simulation"
    if not outputDir.exists():
        outputDir.mkdir()

    # Fatras simulation of muons
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=100000, numThreads=1, logLevel=acts.logging.INFO)

    vertexGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0),
        mean=acts.Vector4(0, 9, 0, 0)
    )

    addParticleGun(
        s=s,
        etaConfig=EtaConfig(10.0, 10.0),
        phiConfig=PhiConfig(0, 0),
        momentumConfig=MomentumConfig(0.5 * u.GeV, 10 * u.GeV),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        vtxGen = vertexGen,
        outputDirRoot=outputDir / "fatras",
    )

    # Simulation algorithm
    alg = acts.examples.FatrasSimulation(
        **acts.examples.defaultKWArgs(
            level=acts.logging.INFO,
            inputParticles="particles_input",
            outputParticlesInitial="particles_initial",
            outputParticlesFinal="particles_final",
            outputSimHits="simhits",
            randomNumbers=rnd,
            trackingGeometry=trackingGeometry,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=True,
            emEnergyLossIonisation=True,
            emEnergyLossRadiation=True,
            emPhotonConversion=True,
            pMin=None,
        )
    )

    # Sequencer
    s.addAlgorithm(alg)

    particlesInitial = alg.config.outputParticlesInitial

    # Only add alias for 'particles_initial' as this is the one we use most
    s.addWhiteboardAlias("particles", particlesInitial)

    s.addAlgorithm(trackEstAlg)

    s.run()

def testLookup():
    # Initialize the geometry
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[4, 10],
        positions=[30, 60, 90],
        stereos=[0, 0, 0],
        binValue=2,
        surfaceType=0
    )

    # Set up the track lookup grid reader
    jsonReaderConfig = acts.examples.JsonTrackLookupGridReader.Config(
        bins=[1, 100],
        xBounds=[-1,1],
        yBounds=[-10,10])
    jsonReader = acts.examples.JsonTrackLookupGridReader(jsonReaderConfig)

    lookupConfig = acts.examples.TrackLookupGridProvider.Config(
        jsonReader,
        "/home/romanurmanov/tools/acts/acts_telescope_geo/ActsTelescopeGeometryDevelopment_build/test-IP.json",
        "/home/romanurmanov/tools/acts/acts_telescope_geo/ActsTelescopeGeometryDevelopment_build/test-REF.json")
    lookup = acts.examples.TrackLookupGridProvider(lookupConfig)

    testConfig = acts.examples.DummyGridTest.Config(
        refLayer=list(trackingGeometry.geoIdSurfaceMap().values())[0],
        lookup=lookup,
        inputHits="simhits",
        inputParticles="particles_input")
    testAlg = acts.examples.DummyGridTest(testConfig, acts.logging.INFO)

    # Set up the magnetic field
    field = acts.ConstantBField(acts.Vector3(50 * u.T, 0, 0))

    outputDir = Path.cwd() / "telescope_simulation"
    if not outputDir.exists():
        outputDir.mkdir()

    # Fatras simulation of muons
    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=1000, numThreads=1, logLevel=acts.logging.INFO)

    vertexGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0),
        mean=acts.Vector4(0, 9, 0, 0)
    )

    addParticleGun(
        s=s,
        etaConfig=EtaConfig(10.0, 10.0),
        phiConfig=PhiConfig(0, 0),
        momentumConfig=MomentumConfig(0.5 * u.GeV, 10 * u.GeV),
        particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        vtxGen = vertexGen,
        outputDirRoot=outputDir / "fatras",
    )

    # Simulation algorithm
    alg = acts.examples.FatrasSimulation(
        **acts.examples.defaultKWArgs(
            level=acts.logging.INFO,
            inputParticles="particles_input",
            outputParticlesInitial="particles_initial",
            outputParticlesFinal="particles_final",
            outputSimHits="simhits",
            randomNumbers=rnd,
            trackingGeometry=trackingGeometry,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=True,
            emEnergyLossIonisation=True,
            emEnergyLossRadiation=True,
            emPhotonConversion=True,
            pMin=None,
        )
    )

    # Sequencer
    s.addAlgorithm(alg)

    particlesInitial = alg.config.outputParticlesInitial

    # Only add alias for 'particles_initial' as this is the one we use most
    s.addWhiteboardAlias("particles", particlesInitial)

    s.addAlgorithm(testAlg)

    s.run()

if "__main__" == __name__:
    # estimateLookup()
    testLookup()