//
//  Connector.swift
//  
//
//  Created by chenlongmingob@gmail.com on 2020/12/28.
//

import Foundation
import CoreBluetooth
import RxSwift


class Connector: NSObject, CBCentralManagerDelegate, CBPeripheralDelegate {
    /// 服务 uuid
    public static let service = "FFE0"
    /// 特征 uuid
    public static let characteristic = "FFE1"
    
    var statePublish = PublishSubject<CMState>()
    public var connecting: CBPeripheral?
    public var connected: CBPeripheral?
    var manager: CBCentralManager!
    
    var service: CBService?
    var characteristic: CBCharacteristic?
    
    public var isConnected: Bool { connected != nil && characteristic != nil }
    
    init(queue: DispatchQueue?, options: [String: Any]?) {
        super.init()
        self.manager = CBCentralManager(delegate: self, queue: queue, options: options)
    }
    
    func startScan(options: [String: Any]? = nil) {
        manager.scanForPeripherals(withServices: nil, options: options)
        statePublish.onNext(.init(state: .startScan))
    }
    
    func stopScan() {
        manager.stopScan()
        statePublish.onNext(.init(state: .stopScan))
    }
    
    func centralManagerDidUpdateState(_ central: CBCentralManager) {
        print("state update: \(central.state)")
    }
    
    func connect(peripheral: CBPeripheral, options: [String: Any]? = nil) {
        connecting = peripheral
        statePublish.onNext(.init(state: .connecting, peripheral: peripheral))
        manager.connect(peripheral, options: options)
    }
    
    func disconnect() {
        if let peripheral = connected {
            manager.cancelPeripheralConnection(peripheral)
        }
    }
    
    func writeData(data: Data) {
        guard let peripheral = connected, let characteristic = characteristic else {
            print("DEBUG: Falha ao escrever: dispositivo ou característica não disponível")
            statePublish.onError(CMError.peripheralDisconnect)
            return
        }
        print("DEBUG: Enviando dados para o dispositivo: \(data.hexEncodedString())")
        let count = (Double(data.count) / 20).rounded(.up)
        let dataEnd = data.count
        for c in 0 ..< Int(count) {
            let start = c * 20
            let chunkEnd = (c + 1) * 20
            let end = chunkEnd < dataEnd ? chunkEnd : dataEnd
            let chunk = data.subdata(in: start ..< end)
            print("DEBUG: Enviando chunk \(c+1): \(chunk.hexEncodedString())")
            peripheral.writeValue(chunk, for: characteristic, type: .withoutResponse)
        }
    }

    
    func centralManager(_ central: CBCentralManager, didDiscover peripheral: CBPeripheral, advertisementData: [String : Any], rssi RSSI: NSNumber) {
        if peripheral.name != nil { print("found:", peripheral.name) }
        statePublish.onNext(.init(state: .scanned, peripheral: peripheral))
    }
    
    func centralManager(_ central: CBCentralManager, didConnect peripheral: CBPeripheral) {
        connecting = nil
        connected = peripheral
        peripheral.delegate = self
        connected?.discoverServices([CBUUID(string: Self.service)])
    }
    
    func centralManager(_ central: CBCentralManager, didDisconnectPeripheral peripheral: CBPeripheral, error: Error?) {
        if let connected = connected, connected == peripheral {
            self.connected = nil
            connecting = nil
            service = nil
            characteristic = nil
            statePublish.onNext(.init(state: .disconnect, peripheral: peripheral))
        }
    }
    
    func centralManager(_ central: CBCentralManager, didFailToConnect peripheral: CBPeripheral, error: Error?) {
        if let connected = connected, connected == peripheral {
            self.connected = nil
        }
        
        if let connecting = connecting, connecting == peripheral {
            self.connecting = nil
        }
        
        statePublish.onError(CMError.failToConnect)
    }
    
    func peripheral(_ peripheral: CBPeripheral, didDiscoverCharacteristicsFor service: CBService, error: Error?) {
        if let characteristics = service.characteristics, let characteristic = characteristics.first(where: { $0.uuid.uuidString == Self.characteristic }) {
            self.characteristic = characteristic
            peripheral.setNotifyValue(true, for: characteristic)
            print("DEBUG: Característica \(Self.characteristic) descoberta e notificação ativada")
            statePublish.onNext(.init(state: .connected, peripheral: peripheral))
        } else {
            print("DEBUG: Falha ao descobrir característica \(Self.characteristic)")
            statePublish.onError(CMError.failToConnect)
        }
    }

    
    func peripheral(_ peripheral: CBPeripheral, didDiscoverServices error: Error?) {
        if let services = peripheral.services, let service = services.first(where: { $0.uuid.uuidString == Self.service }) {
            self.service = service
            peripheral.discoverCharacteristics([CBUUID(string: Self.characteristic)], for: service)
        }
    }
    

    func peripheral(_ peripheral: CBPeripheral, didUpdateValueFor characteristic: CBCharacteristic, error: Error?) {
        statePublish.onNext(.init(state: .notification, data: characteristic.value, peripheral: peripheral))
    }
}

extension Data {
    func hexEncodedString() -> String {
        return map { String(format: "%02x", $0) }.joined()
    }
}

extension CBManagerState: CustomStringConvertible {
    public var description: String {
        switch self {
        case .poweredOff:
            return "poweredOff"
        case .unknown:
            return "unknown"
        case .resetting:
            return "resetting"
        case .unsupported:
            return "unsupported"
        case .unauthorized:
            return "unauthorized"
        case .poweredOn:
            return "poweredOn"
        @unknown default:
            return "not matched"
        }
    }
}
