#include <Arduino.h>

// Pines (ajusta según tu conexión)
const int EOS_PIN = A0;   // Eos
const int TRG_PIN = 3;    // Trigger (conectado a interrupción INT0)
const int SI_PIN = A2;    // Start Integration
const int CLK_PIN = A3;   // Clock (generado por Arduino)
const int VIDEO_PIN = A4; // Salida analógica
const int LED = 2;        // Led

// Variables
volatile bool capture_flag = false;
int iterations = 1;
int delayTime = 1; // delay time
volatile int pixel_count = 0;
volatile bool data_ready = false;
int spectrum[288];         // Array para almacenar el espectro

#define SPEC_CHANNELS    288 // New Spec Channel

void setup() {
  pinMode(EOS_PIN,   OUTPUT);
  pinMode(TRG_PIN,   INPUT);
  pinMode(SI_PIN,    OUTPUT);
  pinMode(CLK_PIN,   OUTPUT);
  pinMode(VIDEO_PIN, INPUT);
  pinMode(LED,       OUTPUT);

  digitalWrite(LED, LOW);
  digitalWrite(CLK_PIN, LOW); // Set CLK_PIN Low
  digitalWrite(SI_PIN, LOW); // Set SI_PIN Low

  // Configurar interrupción por flanco de subida de TRG
  attachInterrupt(digitalPinToInterrupt(TRG_PIN), read_pixel, RISING);

  Serial.begin(115200);
}

int iter = 0;

int enable_read_analog = 0;

int mode_debug = 1; // 0: Usamos Ide Arduino para ver puerto serie
                    // 1: Modo normal aplicación Python

void loop() {
if (mode_debug == 1) {
  if (Serial.available()) {
    String cmd = Serial.readStringUntil('\n');
    cmd.trim();
    
    if (cmd.startsWith("SETTING:")) {
      // Parsear delay e iterations (ej: "SETTING:100,10")
      int commaIndex = cmd.indexOf(',');
      delayTime = cmd.substring(8, commaIndex).toInt();  // "SETTING:" tiene 8 caracteres
      iterations = cmd.substring(commaIndex + 1).toInt();
    
      Serial.print("READY:");  // Confirmar recepción
      Serial.print(delayTime);
      Serial.print(",");
      Serial.println(iterations);
      capture_flag = true;
    }
    else if (cmd == "GET DATA" && capture_flag) {
      enable_read_analog = 1;
      
      for (int i = 0; i < iterations; i++) {
         start_measurement();
         while (!data_ready); // Esperar a que se complete la lectura
         send_spectrum();
        
         if (i < iterations-1) {
           delayMicroseconds(delayTime);
           while (!Serial.available() || Serial.readStringUntil('\n') != "NEXT");
         }
      }
      Serial.println("END DATA");
    }
    else if (cmd == "OSCILLOSCOPE") {
      enable_read_analog = 0;
            
      while (1) {
        start_measurement();
        while (!data_ready); // Esperar a que se complete la lectura         
        
        if (Serial.available()) {
          String cmdosc = Serial.readStringUntil('\n');
          cmdosc.trim();
          if (cmdosc == "STOP OSCILLOSCOPE") {
            // Confirmar recepción comando
            Serial.println("STOP OSCILLOSCOPE");
            
            // Salir del modo OSCILLOSCOPE
            break;
          }
        }
      }  
    }
  }
}
else {
  while (iter < 10) {
    start_measurement();
    if (iter > 0) {
      while (!data_ready); // Esperar a que se complete la lectura
  
      //send_spectrum();
    }
    iter++;
    if (iter == 10)
       Serial.println("END DATA");
  }   
}   
}

void start_measurement() {
  // Resetear variables
  pixel_count = 0;
  data_ready = false;

  // Start clock cycle and set start pulse to signal start
  digitalWrite(CLK_PIN, LOW);
  delayMicroseconds(delayTime);
  digitalWrite(CLK_PIN, HIGH);
  delayMicroseconds(delayTime);
  digitalWrite(CLK_PIN, LOW);
  digitalWrite(SI_PIN, HIGH);
  delayMicroseconds(delayTime);

  //Sample for a period of time
  for(int i = 0; i < 15; i++){

      digitalWrite(CLK_PIN, HIGH);
      delayMicroseconds(delayTime);
      digitalWrite(CLK_PIN, LOW);
      delayMicroseconds(delayTime); 
 
  }

  //Set SI_PIN to low
  digitalWrite(SI_PIN, LOW);
  //Serial.println("START INIT");

  //Sample for a period of time
  for(int i = 0; i < 85; i++){

      digitalWrite(CLK_PIN, HIGH);
      delayMicroseconds(delayTime);
      digitalWrite(CLK_PIN, LOW);
      delayMicroseconds(delayTime); 
      
  }

  //One more clock pulse before the actual read
  digitalWrite(CLK_PIN, HIGH);
  delayMicroseconds(delayTime);
  digitalWrite(CLK_PIN, LOW);
  delayMicroseconds(delayTime);

  //Read from SPEC_VIDEO
  for(int i = 0; i < SPEC_CHANNELS; i++){
      
      digitalWrite(CLK_PIN, HIGH);
      delayMicroseconds(delayTime);      
      digitalWrite(CLK_PIN, LOW);
      delayMicroseconds(delayTime);
        
  }

  //Set SI_PIN to high
  digitalWrite(SI_PIN, HIGH);

  //Sample for a small amount of time
  for(int i = 0; i < 7; i++){
    
      digitalWrite(CLK_PIN, HIGH);
      delayMicroseconds(delayTime);
      digitalWrite(CLK_PIN, LOW);
      delayMicroseconds(delayTime);
    
  }

  digitalWrite(CLK_PIN, HIGH);
  delayMicroseconds(delayTime);
}

// Interrupción: Lee VIDEO_PIN en cada flanco de TRG (a partir del 89º pulso)
void read_pixel() {  
  //Serial.print(analogRead(VIDEO_PIN));
  //Serial.print(",");
  //Serial.println(pixel_count);
  if (pixel_count >= 104 && pixel_count < 104 + 288) {
    int pixel_index = pixel_count - 104;
    if (enable_read_analog == 1)
      spectrum[pixel_index] = analogRead(VIDEO_PIN);
  }
  pixel_count++;

  if (pixel_count == 104 + 288) {
    data_ready = true;  // Lectura completa
  }
}

void send_spectrum() {
  Serial.print("DATA:");
  for (int i = 0; i < 288; i++) {
    Serial.print(spectrum[i]);
    if (i < 287) Serial.print(",");
  }
  Serial.println();
}
