import atmosphere
from constants import *

class Stage:
    def __init__(self, wetMass, dryMass, burnTime, seaISP, vacuumISP):
        self.wetMass = wetMass
        self.dryMass = dryMass
        self.burnTime = burnTime
        self.seaISP = seaISP
        self.vacuumISP = vacuumISP

    @property
    def fuelMass(self):
        return self.wetMass - self.dryMass

    @property
    def massFlowRate(self):
        return self.fuelMass / self.burnTime if self.burnTime > 0 else 0

    @classmethod
    def from_dict(cls, data):
        return cls(wetMass=data["wetMass"], dryMass=data["dryMass"], burnTime=data["burnTime"], seaISP=data["seaISP"], vacuumISP=data["vacuumISP"])

class Rocket:
    def __init__(self, name, stages, dragCoefficient, frontalArea, launchAzimuth, VERTICAL_CLIMB_TIME, PITCH_KICK_DEGREES, PITCH_KICK_DURATION):
        self.name = name
        self.stages = stages
        self.dragCoefficient = dragCoefficient
        self.frontalArea = frontalArea
        self.launchAzimuth = launchAzimuth
        self.VERTICAL_CLIMB_TIME = VERTICAL_CLIMB_TIME
        self.PITCH_KICK_DEGREES = PITCH_KICK_DEGREES
        self.PITCH_KICK_DURATION = PITCH_KICK_DURATION

        self.currentStageIndex = 0
        self.stageTime = 0

        self.mass = sum(stage.wetMass for stage in self.stages)

    @property
    def currentStage(self):
        return self.stages[self.currentStageIndex]

    def calculateThrust(self, height):
        stage = self.currentStage
        pressure_ratio = atmosphere.density(height) / 1.225
        isp = stage.seaISP + (stage.vacuumISP - stage.seaISP) * (1 - pressure_ratio)
        return stage.massFlowRate * isp * SEA_LEVEL_GRAVITY

    def next_stage(self):
        if self.currentStageIndex + 1 < len(self.stages):
            print(f"Stage {self.currentStageIndex + 1} seperation")
            self.mass -= self.currentStage.dryMass
            self.currentStageIndex += 1
            self.stageTime = 0
            return True
        return False

    def update(self, dt, height):
        self.stageTime += dt
        stage = self.currentStage

        if self.stageTime >= stage.burnTime:
            if not self.next_stage():
                return 0
        
        self.mass -= self.currentStage.massFlowRate * dt
        return self.calculateThrust(height)
            
    @classmethod
    def from_dict(cls, data):
        stages = [Stage.from_dict(s) for s in data["stages"]]
        return cls(name=data["name"], stages=stages, dragCoefficient=data["dragCoefficient"], frontalArea=data["frontalArea"], launchAzimuth=data["launchAzimuth"], VERTICAL_CLIMB_TIME=data["VERTICAL_CLIMB_TIME"], PITCH_KICK_DEGREES=data["PITCH_KICK_DEGREES"], PITCH_KICK_DURATION=data["PITCH_KICK_DURATION"])