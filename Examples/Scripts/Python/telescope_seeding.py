#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
)

u = acts.UnitConstants

if "__main__" == __name__:
    detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
        bounds=[200, 200],
        positions=[30, 60, 90, 120, 150, 180, 210, 240, 270],
        stereos=[0, 0, 0, 0, 0, 0, 0, 0, 0],
        binValue=2,
    )

    trackEstConfig = acts.examples.SimTrackParamsEstimationAlgorithm.Config(
        refLayerIds=[list(trackingGeometry.geoIdSurfaceMap().keys())[0]],
        inputHits="simhits",
        inputParticles="particles_input",
        outputIPTrackParameters = "IPTrackParameters",
        outputRefLayerTrackParameters = "RefTrackParameters"
    )

    trackEstAlg = acts.examples.SimTrackParamsEstimationAlgorithm(trackEstConfig, acts.logging.INFO)

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    outputDir = Path.cwd() / "telescope_simulation"
    if not outputDir.exists():
        outputDir.mkdir()

    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(events=1, numThreads=1, logLevel=acts.logging.INFO)

    addParticleGun(
        s,
        EtaConfig(-10.0, 10.0),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(1000, acts.PdgParticle.eMuon, False),
        multiplicity=1,
        rnd=rnd,
        outputDirRoot=outputDir / "fatras",
    )

    # Simulation
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
