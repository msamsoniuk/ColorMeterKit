//
//  CMColor.swift
//  
//
//  Created by chenlongmingob@gmail.com on 2020/12/31.
//

import Foundation

public class CMColor {
    public typealias XYZ = (Double, Double, Double)
    public typealias Lab = (Double, Double, Double)
    public typealias CH = (Double, Double)
    public typealias RGB = (UInt8, UInt8, UInt8)
    
    public var filledSpectral: [Double]!
    private var lightSource: LightSource!
    
    public var xyz: XYZ!
    public var lab: Lab!
    public var ch: CH!
    public var rgb: RGB!
    
    
    public init(spectral: [Double], waveStart: Int, lightSource: LightSource) {
        self.lightSource = lightSource
        fillSpectral(spectral: spectral, waveStart: waveStart)
        let xyz = Self.spectral2XYZ(spectral: filledSpectral, lightSource: lightSource)
        let lab = Self.XYZ2Lab(xyz: xyz, lightSource: lightSource)
        let ch = Self.Lab2CH(lab: lab)
        let rgb = Self.Lab2RGB(lab: lab)
        self.xyz = xyz
        self.lab = lab
        self.ch = ch
        self.rgb = rgb
    }
    
    
    private func fillSpectral(spectral: [Double], waveStart: Int) {
        var spectral = spectral
        spectral.insert(contentsOf: Array(repeating: 0, count: (waveStart - 360) / 10), at: 0)
        spectral.append(contentsOf: Array(repeating: 0, count: (43 - spectral.count)))
        filledSpectral = spectral
    }
    
    
    public static func spectral2XYZ(spectral: [Double], lightSource: LightSource) -> XYZ {
        let params = lightSource.category.kl
        var x: Double = 0
        var y: Double = 0
        var z: Double = 0
        let k = lightSource.kal
        
        for i in 0 ..< 43 {
            let light = params.count <= i ? 0 : params[i]
            let kx = lightSource.angle.kx.count <= i ? 0 : lightSource.angle.kx[i]
            let ky = lightSource.angle.ky.count <= i ? 0 : lightSource.angle.ky[i]
            let kz = lightSource.angle.kz.count <= i ? 0 : lightSource.angle.kz[i]
            
            x += spectral[i] * 0.01 * light * kx;
            y += spectral[i] * 0.01 * light * ky;
            z += spectral[i] * 0.01 * light * kz;
        }
    
        return (x * k, y * k, z * k)
    }
    
    
    public static func XYZ2Lab(xyz: XYZ, lightSource: LightSource) -> Lab {
        var l: Double = 0
        var a: Double = 0
        var b: Double = 0
        let (x, y, z) = xyz
        
        var paramX: Double = x / lightSource.kxyz_labch[0]
        var paramY: Double = y / lightSource.kxyz_labch[1]
        var paramZ: Double = z / lightSource.kxyz_labch[2]
        
        if paramX > 0.008856 {
           paramX = pow(paramX, 0.3333333)
        } else {
            paramX = 7.787 * paramX + 0.138
        }
        
        if (paramY > 0.008856) {
            paramY = pow(paramY, 0.3333333);
            l = 116 * paramY - 16;
        } else {
            l = 903.3 * paramY;
            paramY = (7.787 * paramY) + 0.138;
        }
        
        if paramZ > 0.008856 {
           paramZ = pow(paramZ, 0.3333333)
        } else {
            paramZ = 7.787 * paramZ + 0.138
        }
        
        a = 500 * (paramX - paramY)
        b = 200 * (paramY - paramZ)
        
        return (l, a, b)
    }
    
    
    public static func Lab2CH(lab: Lab) -> CH {
        var (l, a, b) = lab
        if l < 0 { l = 0 }
        let c = sqrt(pow(a, 2) + pow(b, 2))
        var h = 0.0
        if a == 0 && b > 0 {
            h = 90
        } else if a == 0 && b < 0 {
            h = 270
        } else if a >= 0 && b == 0 {
            h = 0
        } else if a < 0 && b == 0 {
            h = 180
        } else {
            h = atan(b / a) * 57.3
            if a > 0 && b > 0 {  }
            else if a < 0 {
                h += 180
            } else {
                h += 360
            }
        }
        
        return (c, h)
    }
    
    
    public static func Lab2RGB(lab: Lab) -> RGB {
        let (l, a, b) = lab
        
        var y = (l + 16) / 116
        var x = a / 500 + y
        var z = y - b / 200
        
        y = y > 6 / 29 ? pow(y, 3) : (y - 16 / 116) / 7.787
        x = x > 6 / 29 ? pow(x, 3) : (x - 16 / 116) / 7.787
        z = z > 6 / 29 ? pow(z, 3) : (z - 16 / 116) / 7.787
        
        x *= 0.95047
        z *= 1.08883
        
        var red = 3.2406 * x - 1.5372 * y - 0.4986 * z
        var green = -0.9689 * x + 1.8758 * y + 0.0415 * z
        var blue = 0.0557 * x - 0.2040 * y + 1.0570 * z
        
        red = red > 0.0031308 ? 1.055 * pow(red, 1 / 2.4) - 0.055 : 12.92 * red
        green = green > 0.0031308 ? 1.055 * pow(green, 1 / 2.4) - 0.055 : 12.92 * green
        blue = blue > 0.0031308 ? 1.055 * pow(blue, 1 / 2.4) - 0.055 : 12.92 * blue
        
        let R = max(min(255, red * 255), 0)
        let G = max(min(255, green * 255), 0)
        let B = max(min(255, blue * 255), 0)
        
        
        return (UInt8(R), UInt8(G), UInt8(B))
    }


    public static func Lab2RGBLinear(lab: (Double, Double, Double)) -> (Double, Double, Double) {
        let (L, a, b) = lab

        // XYZ reference white (D65)
        let Xn = 0.95047
        let Yn = 1.00000
        let Zn = 1.08883

        // Convert LAB to XYZ
        let fy = (L + 16.0) / 116.0
        let fx = a / 500.0 + fy
        let fz = fy - b / 200.0

        func fInv(_ t: Double) -> Double {
            let cube = t * t * t
            return cube > 0.008856 ? cube : (t - 16.0 / 116.0) / 7.787
        }

        let x = Xn * fInv(fx)
        let y = Yn * fInv(fy)
        let z = Zn * fInv(fz)

        // Convert XYZ to linear RGB
        var rl =  3.2406 * x - 1.5372 * y - 0.4986 * z
        var gl = -0.9689 * x + 1.8758 * y + 0.0415 * z
        var bl =  0.0557 * x - 0.2040 * y + 1.0570 * z

        // Clamp to 0.0–1.0 range before applying gamma
        rl = max(0.0, min(1.0, rl))
        gl = max(0.0, min(1.0, gl))
        bl = max(0.0, min(1.0, bl))

        // Apply gamma correction (sRGB)
        func gammaCorrect(_ c: Double) -> Double {
            return c <= 0.0031308 ? 12.92 * c : 1.055 * pow(c, 1/2.4) - 0.055
        }

        let rr = gammaCorrect(rl)
        let gg = gammaCorrect(gl)
        let bb = gammaCorrect(bl)

        return (rr, gg, bb) // RGB values in 0.0–1.0 range
    }

    public static func rgbToLab(_ rgb: (UInt8, UInt8, UInt8)) -> (Double, Double, Double) {
            let R = Double(rgb.0) / 255.0
            let G = Double(rgb.1) / 255.0
            let B = Double(rgb.2) / 255.0

            // Linearize
            func invGamma(_ c: Double) -> Double {
                return c <= 0.04045 ? c / 12.92 : pow((c + 0.055) / 1.055, 2.4)
            }

            let rl = invGamma(R)
            let gl = invGamma(G)
            let bl = invGamma(B)

            // RGB → XYZ (sRGB, D65)
            let x = rl * 0.4124 + gl * 0.3576 + bl * 0.1805
            let y = rl * 0.2126 + gl * 0.7152 + bl * 0.0722
            let z = rl * 0.0193 + gl * 0.1192 + bl * 0.9505

            // Normalize with D65 white point
            let Xn = 0.95047
            let Yn = 1.00000
            let Zn = 1.08883

            var fx = x / Xn
            var fy = y / Yn
            var fz = z / Zn

            func f(_ t: Double) -> Double {
                return t > 0.008856 ? pow(t, 1.0/3.0) : (7.787 * t) + (16.0 / 116.0)
            }

            fx = f(fx)
            fy = f(fy)
            fz = f(fz)

            let L = (116 * fy) - 16
            let a = 500 * (fx - fy)
            let b = 200 * (fy - fz)

            return (L, a, b)
        }

}

extension CMColor {
    public typealias CMYK = (Double, Double, Double, Double)
    
    public var cmyk: CMYK {
        return Self.RGB2CMYK(rgb: rgb)
    }
    
    public static func RGB2CMYK(rgb: RGB) -> CMYK {
        // Normaliza RGB de 0-255 para 0-1
        let r = Double(rgb.0) / 255.0
        let g = Double(rgb.1) / 255.0
        let b = Double(rgb.2) / 255.0
        
        // Calcula K (preto)
        let k = 1.0 - max(r, g, b)
        
        // Evita divisão por zero
        let divisor = (1.0 - k) == 0 ? 1.0 : (1.0 - k)
        
        // Calcula C, M, Y
        let c = (1.0 - r - k) / divisor
        let m = (1.0 - g - k) / divisor
        let y = (1.0 - b - k) / divisor
        
        // Converte para porcentagem (0-100)
        return (c * 100.0, m * 100.0, y * 100.0, k * 100.0)
    }
}
