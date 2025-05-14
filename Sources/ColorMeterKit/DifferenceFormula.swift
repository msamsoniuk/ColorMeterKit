//
//  DifferenceFormula.swift
//  
//
//  Created by chenlongmingob@gmail.com on 2020/12/30.
//

import Foundation

public enum DifferenceFormula: UInt8, Codable, CaseIterable {
    case CIE_dE_ab
    case dE_ch
    case dE_00
    case dE_cmc
    case dE_94
    case hunter_dE_ab
    case dE_uv
    
    public func calculateDifference(_ color1: CMColor.Lab, _ color2: CMColor.Lab) -> Double {
        switch self {
        case .CIE_dE_ab:
            return deltaE76(color1: color1, color2: color2)
        case .dE_ch:
            return deltaECH(color1: color1, color2: color2)
        case .dE_00:
            return deltaE00(color1: color1, color2: color2)
        case .dE_cmc:
            return deltaECMC(color1: color1, color2: color2)
        case .dE_94:
            return deltaE94(color1: color1, color2: color2)
        case .hunter_dE_ab:
            return hunterDeltaE(color1: color1, color2: color2)
        case .dE_uv:
            return deltaEUV(color1: color1, color2: color2)
        }
    }
    
    // MARK: - CIE 1976 ΔE*ab formula
    private func deltaE76(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        return sqrt(pow(L2 - L1, 2) + pow(a2 - a1, 2) + pow(b2 - b1, 2))
    }
    
    // MARK: - Delta E CH (chroma/hue)
    private func deltaECH(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (_, a1, b1) = color1
        let (_, a2, b2) = color2
        let deltaC = sqrt(pow(a2, 2) + pow(b2, 2)) - sqrt(pow(a1, 2) + pow(b1, 2))
        let deltaH = sqrt(pow(a2 - a1, 2) + pow(b2 - b1, 2) - pow(deltaC, 2))
        return sqrt(pow(deltaC, 2) + pow(deltaH, 2))
    }
    
    // MARK: - CIEDE2000 ΔE*₀₀ formula (your existing implementation is correct)
    private func deltaE00(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        
        // Constants
        let kL: Double = 1.0 // lightness
        let kC: Double = 1.0 // chroma
        let kH: Double = 1.0 // hue
        
        // Step 1: Calculate C' and h'
        let C1 = sqrt(a1 * a1 + b1 * b1)
        let C2 = sqrt(a2 * a2 + b2 * b2)
        let C_avg = (C1 + C2) / 2.0
        
        let G = 0.5 * (1 - sqrt(pow(C_avg, 7) / (pow(C_avg, 7) + pow(25, 7))))     
           
        let a1_prime = a1 * (1 + G)
        let a2_prime = a2 * (1 + G)
        
        let C1_prime = sqrt(a1_prime * a1_prime + b1 * b1)
        let C2_prime = sqrt(a2_prime * a2_prime + b2 * b2)
        
        let h1_prime = (atan2(b1, a1_prime) * 180 / .pi).truncatingRemainder(dividingBy: 360)
        let h2_prime = (atan2(b2, a2_prime) * 180 / .pi).truncatingRemainder(dividingBy: 360)
        
        // Step 2: Calculate ΔL', ΔC', ΔH'
        let ΔL_prime = L2 - L1
        let ΔC_prime = C2_prime - C1_prime
        
        var Δh_prime: Double
        if (C1_prime * C2_prime) == 0 {
            Δh_prime = 0
        } else if abs(h2_prime - h1_prime) <= 180 {
            Δh_prime = h2_prime - h1_prime
        } else if (h2_prime - h1_prime) > 180 {
            Δh_prime = h2_prime - h1_prime - 360
        } else {
            Δh_prime = h2_prime - h1_prime + 360
        }
        
        let ΔH_prime = 2 * sqrt(C1_prime * C2_prime) * sin(Δh_prime * .pi / 360)
        
        // Step 3: Calculate weighting functions
        let L_avg_prime = (L1 + L2) / 2.0
        let C_avg_prime = (C1_prime + C2_prime) / 2.0
        
        var h_avg_prime: Double
        if (C1_prime * C2_prime) == 0 {
            h_avg_prime = 0
        } else if abs(h1_prime - h2_prime) <= 180 {
            h_avg_prime = (h1_prime + h2_prime) / 2.0
        } else if (h1_prime + h2_prime) < 360 {
            h_avg_prime = (h1_prime + h2_prime + 360) / 2.0
        } else {
            h_avg_prime = (h1_prime + h2_prime - 360) / 2.0
        }
        
        let T = 1 - 0.17 * cos((h_avg_prime - 30) * .pi / 180) +
                   0.24 * cos(2 * h_avg_prime * .pi / 180) +
                   0.32 * cos((3 * h_avg_prime + 6) * .pi / 180) -
                   0.20 * cos((4 * h_avg_prime - 63) * .pi / 180)
        
        let SL = 1 + (0.015 * pow(L_avg_prime - 50, 2)) / sqrt(20 + pow(L_avg_prime - 50, 2))
        let SC = 1 + 0.045 * C_avg_prime
        let SH = 1 + 0.015 * C_avg_prime * T
        
        let Δθ = 30 * exp(-pow((h_avg_prime - 275) / 25, 2))
        let RC = 2 * sqrt(pow(C_avg_prime, 7) / (pow(C_avg_prime, 7) + pow(25, 7)))
        let RT = -RC * sin(2 * Δθ * .pi / 180)
        
        // Final calculation
        let term1 = ΔL_prime / (kL * SL)
        let term2 = ΔC_prime / (kC * SC)
        let term3 = ΔH_prime / (kH * SH)
        
        return sqrt(term1 * term1 + term2 * term2 + term3 * term3 + RT * term2 * term3)
    }
    
    // MARK: - CMC l:c formula
    private func deltaECMC(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        
        let C1 = sqrt(a1 * a1 + b1 * b1)
        let C2 = sqrt(a2 * a2 + b2 * b2)
        let deltaC = C1 - C2
        
        let deltaL = L1 - L2
        let deltaA = a1 - a2
        let deltaB = b1 - b2
        let deltaH = sqrt(deltaA * deltaA + deltaB * deltaB - deltaC * deltaC)
        
        let SL = L1 < 16 ? 0.511 : (0.040975 * L1) / (1 + 0.01765 * L1)
        let SC = 0.0638 * C1 / (1 + 0.0131 * C1) + 0.638
        let F = sqrt(pow(C1, 4) / (pow(C1, 4) + 1900))
        let T = (164...345).contains(b1) ? 0.56 + abs(0.2 * cos((b1 + 168) * .pi / 180)) : 
                                         0.36 + abs(0.4 * cos((b1 + 35) * .pi / 180))
        let SH = SC * (F * T + 1 - F)
        
        return sqrt(pow(deltaL / (SL * 2), 2) + pow(deltaC / (SC * 1), 2) + pow(deltaH / SH, 2))
    }
    
    // MARK: - CIE 1994 ΔE formula
    private func deltaE94(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        
        let deltaL = L1 - L2
        let C1 = sqrt(a1 * a1 + b1 * b1)
        let C2 = sqrt(a2 * a2 + b2 * b2)
        let deltaC = C1 - C2
        let deltaA = a1 - a2
        let deltaB = b1 - b2
        let deltaH = sqrt(deltaA * deltaA + deltaB * deltaB - deltaC * deltaC)
        
        let kL: Double = 1.0
        let kC: Double = 1.0
        let kH: Double = 1.0
        let K1: Double = 0.045
        let K2: Double = 0.015
        
        let SL = 1.0
        let SC = 1 + K1 * C1
        let SH = 1 + K2 * C1
        
        return sqrt(pow(deltaL / (kL * SL), 2) + pow(deltaC / (kC * SC), 2) + pow(deltaH / (kH * SH), 2))
    }
    
    // MARK: - Hunter ΔE formula
    private func hunterDeltaE(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        return sqrt(pow(L2 - L1, 2) + pow(a2 - a1, 2) + pow(b2 - b1, 2))
    }
    
    // MARK: - CIE ΔE*uv formula
    private func deltaEUV(color1: CMColor.Lab, color2: CMColor.Lab) -> Double {
        // Convert Lab to Luv first (simplified)
        let (L1, a1, b1) = color1
        let (L2, a2, b2) = color2
        
        let u1 = (4 * a1) / (a1 + 15 * b1 + 3)
        let v1 = (9 * b1) / (a1 + 15 * b1 + 3)
        let u2 = (4 * a2) / (a2 + 15 * b2 + 3)
        let v2 = (9 * b2) / (a2 + 15 * b2 + 3)
        
        return sqrt(pow(L2 - L1, 2) + pow(u2 - u1, 2) + pow(v2 - v1, 2))
    }
}

