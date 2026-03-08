import json

def exportTrajectory(position, velocity, planet, filename="trajectory.json"):
    pos = position[::10].tolist()
    vel = velocity[::10].tolist()
    
    data = {
        "planetRadius": planet.RADIUS,
        "position": pos,
        "velocity": vel
    }
    
    with open(filename, "w") as f:
        json.dump(data, f)
    
    print(f"Trajectory exported to {filename}")