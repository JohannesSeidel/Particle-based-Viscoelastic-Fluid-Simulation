#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 19:42:37 2020

Based on the paper: Particle-based Viscoelastic Fluid Simulation
by:                 Simon Clavet, Philippe Beaudoin, and Pierre Poulin


@author: johannesseidel
"""

import numpy as np




# ------------------- Classes & Functions ------------------ #

class Spring:
    def __init__(self, a, b, c=0.0):
        self.a           = a
        self.b           = b
        self.rest_length = c


class Particle:
    def __init__(self, position, velocity):
        self.prev_position = position
        self.postition     = position
        self.velocity      = velocity
        self.spring        = []
    
    
class SPH:
    def __init__(self, end_time=100.0, delta_t=1.0):
        self.gravity  = 9.81
        self.rest_rho = 10.0
        
        self.k         = 0.004
        self.near_k    = 0.01
        self.spring_k  = 0.3
        
        self.alpha     = 0.3
        self.beta      = 1.0
        self.gamma     = 1.0
        self.sigma     = 1.0
        
        self.h         = 1.0 # interaction range
        
        self.delta_t   = delta_t
        self.end_time  = end_time
        
        self.particles = []
        self.springs   = []
    
    
    def apply_gravity(self):
        for particle in self.particles:
            particle.velocity += self.delta_t * self.gravity
    
    
    def apply_viscosity(self):
        for particle_i in self.particles:
            for particle_j in self.particles:
                if particle_i is particle_j:
                    break
                r_ij = np.linalg.norm(particle_i.position - particle_j.position)
                if r_ij < self.h:
                    vector = (particle_i.position - particle_j.position) / r_ij
                    # inward radial velocity
                    u = (particle_i.velocity - particle_j.velocity) * vector
                    if u > 0:
                        # linear and quadratic impulses
                        I = self.delta_t * (1 - r_ij/self.h) * (self.sigma * u + self.beta * u**2) * r_ij
                        particle_i.velocity -= I/2
                        particle_j.velovity += I/2
                        
                        
    def advance_to_predicted_position(self):
        for particle in self.particles:
            # save previous position
            particle.prev_position = particle.position
            # advance to predicted position
            particle.position += self.delta_t * particle.velocity
        
        
    def adjust_springs(self):
        for particle_i in self.particles:
            for particle_j in self.particles:
                if particle_j is particle_i:
                    break
                r_ij = np.linalg.norm(particle_i.position - particle_j.position)
                if r_ij < self.h:
                    idx = [count for count, elem in enumerate(self.springs) if elem.a is particle_i and elem.b is particle_j]
                    if not idx:
                        self.springs.append(Spring(particle_i, particle_j, self.h))
                        idx = -1
                    # tolerable deformation = yield ratio * rest length
                    d = self.gamma * self.springs[idx].rest_length
                    if r_ij > self.springs[idx].rest_length + d: # stretch
                        self.springs[idx].rest_length += self.delta_t * self.alpha * (r_ij - self.springs[idx].rest_length - d)
                    elif r_ij < self.springs[idx].rest_length - d: # compress
                        self.springs[idx].rest_length -= self.delta_t * self.alpha * (self.springs[idx].rest_length - d - r_ij)
        self.springs = [spring for spring in self.springs if spring.rest_length > self.h]
                    
    
    def double_density_relaxation(self):
        for particle in self.particles:
            rho      = 0
            near_rho = 0
            # compute density and near-density
            for neighbour in (x for x in self.particles if x is not particle):
                q = np.linalg.norm(particle.position - neighbour.position) / self.h
                if q < 1:
                    rho      += (1-q)**2
                    near_rho += (1-q)**3
            # compute pressure and near-pressure
            p      = self.k      * (rho - self.rest_rho)
            near_p = self.near_k * near_rho
            dx = np.array([0, 0, 0])
            for neighbour in (x for x in self.particles if x is not particle):
                vector = particle.position - neighbour.position
                q = np.linalg.norm(vector) / self.h
                if q < 1:
                    # apply displacements
                    D = self.delta_t**2 * (p * (1 - q) + near_p * (1 - q)**2) * vector / np.linalg.norm(vector)
                    neighbour.position += D/2
                    dx -= D/2
            particle.position += dx
         
            
    def compute_velocity(self):
        for particle in self.particles:
            # use previous position to compute next velocity
            particle.velocity = (particle.position - particle.prev_position) / self.delta_t    
            
            
    def run(self):
        for t in range(self.end_time/self.delta_t):
            self.apply_gravity()
            self.apply_viscosity()
            self.advance_to_predicted_position()
            self.adjust_springs()
            self.apply_spring_displacements()
            self.double_density_relaxation()
            self.resolve_collisions()
            self.compute_velocity()
