// ===========================
// Optogenetic TTL Stimulation (Serial-Configurable)
// ===========================
// Features
// - Control TTL pulses via serial commands like:
//     20Hz,5ms,10s
//     1 Hz, 10 ms, 1 s
//     40, 2, 5          (positional defaults: Hz, ms, s)
// - Type 'help' for usage. Type 'stop' (or send 'x') during stimulation to abort.
// - Echoes parsed parameters and performs basic validation.
// - Non-blocking serial reads during stimulation to allow abort.
//
// Wiring
//   laser/LED TTL input  <-- pin 5 (laserPin)
//   laser/LED ground     <-- GND
//   (Use a proper laser driver. For LED test on 5V, include ~220 Ω series resistor.)

const uint8_t  laserPin = 5;   // DIGITAL OUTPUT PIN FOR TTL

// "Last used" defaults (used if user just sends 's')
float lastFreqHz = 20.0f;
float lastPulseMs = 5.0f;
float lastDurMs = 10000.0f;

// ----------------- Utility helpers -----------------
void toLowerInPlace(String &s) {
  for (size_t i = 0; i < s.length(); ++i) s[i] = (char)tolower(s[i]);
}

// Extract leading number from token (supports float); returns NAN on failure
float parseLeadingNumber(const String &tok) {
  int i = 0;
  bool seenDigit = false, seenDot = false;
  while (i < tok.length()) {
    char c = tok[i];
    if (isdigit(c)) { seenDigit = true; i++; }
    else if (c == '.' && !seenDot) { seenDot = true; i++; }
    else break;
  }
  if (!seenDigit) return NAN;
  return tok.substring(0, i).toFloat();
}

// Returns true if token contains substring (case-insensitive)
bool containsCI(const String &tok, const char *needle) {
  String t = tok; toLowerInPlace(t);
  return t.indexOf(needle) >= 0;
}

// Split a line by commas or whitespace into up to maxParts tokens
int splitFlexible(String s, String parts[], int maxParts) {
  // Replace commas and semicolons with spaces to unify
  for (int i = 0; i < s.length(); ++i) {
    char c = s[i];
    if (c == ',' || c == ';' || c == '\t') s[i] = ' ';
  }
  s.trim();
  int count = 0;
  int start = 0;
  while (count < maxParts) {
    // Skip spaces
    while (start < s.length() && s[start] == ' ') start++;
    if (start >= s.length()) break;
    int end = start;
    while (end < s.length() && s[end] != ' ') end++;
    parts[count++] = s.substring(start, end);
    start = end + 1;
  }
  return count;
}

// Parse a command line into freq (Hz), pulse (ms), dur (ms)
// Accepted forms:
//  - "20Hz 5ms 10s" (any separators)
//  - "1hz 10 ms 1 s"
//  - "40 2 5" (positional defaults: Hz, ms, s)
// Returns true on success.
bool parseCommandLine(String line, float &freqHz, float &pulseMs, float &durMs) {
  line.trim();
  toLowerInPlace(line);

  if (line == "s" || line == "start") { // use last params
    freqHz = lastFreqHz; pulseMs = lastPulseMs; durMs = lastDurMs; return true;
  }

  String toks[6];
  int n = splitFlexible(line, toks, 6);
  if (n == 0) return false;

  // If user typed a keyword command, handle it elsewhere
  if (toks[0] == "help" || toks[0] == "stop" || toks[0] == "x") return false;

  // Try to assign by units first
  float f = NAN, pw_ms = NAN, dur_ms = NAN;
  for (int i = 0; i < n; ++i) {
    String t = toks[i];
    float val = parseLeadingNumber(t);
    if (isnan(val)) continue;

    if (containsCI(t, "hz")) {
      f = val;
    } else if (containsCI(t, "ms")) {
      // token declares milliseconds explicitly
      if (isnan(pw_ms)) pw_ms = val; else if (isnan(dur_ms)) dur_ms = val; else {/*extra ignored*/}
    } else if (containsCI(t, "sec") || (t.endsWith("s") && !t.endsWith("ms"))) {
      // seconds (e.g., 10s, 10sec, 10seconds)
      float ms = val * 1000.0f;
      if (isnan(dur_ms)) dur_ms = ms; else if (isnan(pw_ms)) pw_ms = ms; else {/*extra ignored*/}
    } else {
      // No unit: assign positionally later
    }
  }

  // Fill in remaining from positional defaults: Hz, ms, s
  int pos = 0;
  for (int i = 0; i < n; ++i) {
    String t = toks[i];
    float val = parseLeadingNumber(t);
    if (isnan(val)) continue;

    bool hasUnit = containsCI(t, "hz") || containsCI(t, "ms") || containsCI(t, "sec") || (t.endsWith("s") && !t.endsWith("ms"));
    if (hasUnit) continue; // already handled

    if (pos == 0 && isnan(f)) { f = val; pos++; }
    else if (pos <= 1 && isnan(pw_ms)) { pw_ms = val; pos++; }
    else if (pos <= 2 && isnan(dur_ms)) { dur_ms = val * 1000.0f; pos++; }
  }

  // Final sanity check / defaults
  if (isnan(f) || isnan(pw_ms) || isnan(dur_ms)) return false;

  freqHz  = f;
  pulseMs = pw_ms;
  durMs   = dur_ms;
  return true;
}

// Deliver TTL pulses; returns when done or aborted.
bool deliverPulses(float freqHz, float pulseMs, float durMs) {
  if (freqHz <= 0 || pulseMs <= 0 || durMs <= 0) return false;

  float periodMs = 1000.0f / freqHz;
  if (pulseMs >= periodMs) {
    Serial.println(F("! Pulse width >= period; clamping to 90% of period."));
    pulseMs = periodMs * 0.9f;
  }

  unsigned long uPulse = (unsigned long)(pulseMs + 0.5f);
  unsigned long uOff   = (unsigned long)(periodMs - pulseMs + 0.5f);
  unsigned long totalPulses = (unsigned long)((durMs / periodMs) + 0.5f);

  Serial.print(F("Pulses: ")); Serial.print(totalPulses);
  Serial.print(F(" | Period (ms): ")); Serial.print(periodMs, 3);
  Serial.print(F(" | High (ms): ")); Serial.print(uPulse);
  Serial.print(F(" | Low (ms): ")); Serial.println(uOff);

  unsigned long started = millis();
  for (unsigned long i = 0; i < totalPulses; ++i) {
    // Abort if user types 'x' or 'stop'
    if (Serial.available()) {
      String s = Serial.readStringUntil('\n');
      s.trim(); toLowerInPlace(s);
      if (s == "x" || s == "stop" || s == "abort") {
        Serial.println(F("<< Aborted by user."));
        digitalWrite(laserPin, LOW);
        return false;
      }
    }

    digitalWrite(laserPin, HIGH);
    delay(uPulse);
    digitalWrite(laserPin, LOW);
    delay(uOff);
  }

  unsigned long elapsed = millis() - started;
  Serial.print(F(">> Stimulation complete. Elapsed ms: "));
  Serial.println(elapsed);
  return true;
}

// ----------------- Arduino setup/loop -----------------
void setup() {
  pinMode(laserPin, OUTPUT);
  digitalWrite(laserPin, LOW);
  Serial.begin(115200);
  delay(50);
  Serial.println(F("\nOpto TTL ready."));
  Serial.println(F("Type examples: 20Hz,5ms,10s  |  1 Hz 10 ms 1 s  |  40 2 5"));
  Serial.println(F("Commands: 's' to start with last params, 'stop' or 'x' to abort, 'help' for help."));
}

void loop() {
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    line.trim();
    String lower = line; toLowerInPlace(lower);

    if (lower == "help") {
      Serial.println(F("Usage: <freq><Hz>, <pulse><ms>, <duration><s>"));
      Serial.println(F("Examples: 20Hz 5ms 10s | 1hz,10ms,1s | 40 2 5 (defaults Hz,ms,s)"));
      Serial.println(F("Then stimulation starts immediately. Type 'x' to abort."));
      return;
    }

    float f, pw, dur;
    if (!parseCommandLine(line, f, pw, dur)) {
      Serial.println(F("! Could not parse. Type 'help' for examples."));
      return;
    }

    // Save as last used
    lastFreqHz = f; lastPulseMs = pw; lastDurMs = dur;

    Serial.print(F(">> Starting: ")); Serial.print(f, 3); Serial.print(F(" Hz, "));
    Serial.print(pw, 3); Serial.print(F(" ms, "));
    Serial.print(dur/1000.0f, 3); Serial.println(F(" s"));

    deliverPulses(f, pw, dur);
  }
}
