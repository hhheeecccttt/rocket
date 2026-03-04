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
        return cls(
            wetMass=data["wetMass"],
            dryMass=data["dryMass"],
            burnTime=data["burnTime"],
            seaISP=data["seaISP"],
            vacuumISP=data["vacuumISP"]
        )

class Rocket:
    def __init__(self, name, stages, dragCoefficient, frontalArea,
                 VERTICAL_CLIMB_TIME, PITCH_KICK_DEGREES, PITCH_KICK_DURATION):
        
        self.name = name
        self.stages = stages  # list of Stage objects
        self.dragCoefficient = dragCoefficient
        self.frontalArea = frontalArea
        self.VERTICAL_CLIMB_TIME = VERTICAL_CLIMB_TIME
        self.PITCH_KICK_DEGREES = PITCH_KICK_DEGREES
        self.PITCH_KICK_DURATION = PITCH_KICK_DURATION

        self.currentStageIndex = 0
        self.stageTime = 0

    @property
    def wetMass(self):
        return sum(stage.wetMass for stage in self.stages)

    @property
    def currentStage(self):
        return self.stages[self.currentStageIndex]

    def next_stage(self):
        if self.currentStageIndex + 1 < len(self.stages):
            print(f"Stage {self.currentStageIndex + 1} sep")
            self.currentStageIndex += 1
            self.stageTime = 0
            return True
        return False

    @classmethod
    def from_dict(cls, data):
        stages = [Stage.from_dict(s) for s in data["stages"]]
        return cls(
            name=data["name"],
            stages=stages,
            dragCoefficient=data["dragCoefficient"],
            frontalArea=data["frontalArea"],
            VERTICAL_CLIMB_TIME=data["VERTICAL_CLIMB_TIME"],
            PITCH_KICK_DEGREES=data["PITCH_KICK_DEGREES"],
            PITCH_KICK_DURATION=data["PITCH_KICK_DURATION"]
        )