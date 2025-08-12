// ===========================
// Optogenetic TTL Stimulation
// ===========================

const uint8_t laserPin    = 5;      // TTL output pin
const uint16_t frequency  = 20;     // pulse frequency in Hz
const uint16_t pulseWidth = 5;      // pulse width in ms
const uint32_t duration   = 10000;  // total stimulation time in ms

void setup() {
  pinMode(laserPin, OUTPUT);
  digitalWrite(laserPin, LOW);
  Serial.begin(9600);
  Serial.println("Ready. Type 's' (no line ending) and press Enter to start.");
}

void loop() {
  // Check for incoming serial character
  if (Serial.available() > 0) {
    char c = Serial.read();
    if (c == 's') {
      Serial.println(">> Stimulation started");
      
      // Compute period and number of pulses
      const uint16_t period = 1000 / frequency;                   // ms
      const uint16_t offTime = period - pulseWidth;               // ms
      const uint16_t totalPulses = frequency * duration / 1000;   // e.g. 20 Hz * 10 s = 200 pulses

      // Deliver pulses
      for (uint16_t i = 0; i < totalPulses; ++i) {
        digitalWrite(laserPin, HIGH);
        delay(pulseWidth);
        digitalWrite(laserPin, LOW);
        delay(offTime);
      }

      Serial.println(">> Stimulation ended");
    }
  }
}
