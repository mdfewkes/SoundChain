#pragma once

#include "SoundChain.hpp"
#include <vector>
#include <array>
#include <atomic>

namespace TBL {

	struct WaveTable {
		static const int DEFAULT_TABLE_SIZE = 1024;

		std::vector<float> table = {};
		int tableSize = DEFAULT_TABLE_SIZE;

		void BuildTable() {
			table.resize(tableSize, 0.0f);

			for (int position = 0; position < tableSize; position++) {
				double phase = (double)position / (double)tableSize;
				table[position] = sin(phase * 2 * M_PI);
			}
		}

		float Lookup(float phase) {
			phase = (phase - floor(phase)) * tableSize;

			float t = phase - floor(phase);
			int startPosition = floor(phase);
			int endPosition = startPosition + 1;
			if (endPosition >= tableSize) endPosition -= tableSize;

			return (t * table[endPosition]) + ((1.0f - t) * table[startPosition]);
		}

		float LookupSin(float phase) {
			phase = phase - floor(phase);
			return std::sin(M_PI * 2.0 * phase);
		}
	};

	template <typename T, int S>
	class RingBuffer {
	public:
		static_assert(S > 0, "Size must be positive");

		bool Push(const T& item) {
			int headIndex = head.load(std::memory_order_relaxed);
			int next = NextIndex(headIndex);
			if (next == tail.load(std::memory_order_acquire)) {
				return false; // full
			}
			buffer[headIndex] = item;
			head.store(next, std::memory_order_release);
			return true;
		}

		T Pop() {
			T item;

			int tailIndex = tail.load(std::memory_order_relaxed);
			if (tailIndex == head.load(std::memory_order_acquire)) {
				return item; // empty
			}
			item = buffer[tailIndex];
			tail.store(NextIndex(tailIndex), std::memory_order_release);
			return item;
		}

		bool isEmpty() const {
			return head.load(std::memory_order_acquire) ==
				   tail.load(std::memory_order_acquire);
		}

	private:
		constexpr int NextIndex(int i) const noexcept {
			return (i + 1) % S;
		}

		std::array<T, S+1> buffer;
		std::atomic<int> head{0};
		std::atomic<int> tail{0};
	};


	class AdditiveOscSoundChain : public SoundChainBase {
	public:
		AdditiveOscSoundChain(float frequency = 440.0f) : _baseFrequency(frequency) {}
		~AdditiveOscSoundChain() {}

		float GetFrequency() {return _baseFrequency;}
		void SetFrequency(float frequency) {
			int offParams = 1 - currentParamIndex.load();

			params[offParams].frequency = frequency;
			if (params[offParams].oscInc.size() == 0) return;

			// float osc_mag = 0;
			params[offParams].numberOfActiveOsc = 0;
			for (int i = 0; i < NUMBER_OF_OSC; i++) {
				params[offParams].oscInc[i] = frequency * params[offParams].oscHarmonic[i] / ReadSettings().SampleRate;

				if (params[offParams].oscInc[i] < 0.5f) {
					params[offParams].numberOfActiveOsc++;
					// osc_mag += params[offParams].oscAmplitude[i];
				} else {
					break;
				}
			}
			// params[offParams].amp = 1.0f / osc_mag * 0.5f;

			_baseFrequency = frequency;
			currentParamIndex.store(offParams);
		}

	private:
		static const int NUMBER_OF_OSC = 32;

		struct SynthParams {
			int numberOfActiveOsc = 0;
			std::array<float, NUMBER_OF_OSC> oscInc = {};
			std::array<float, NUMBER_OF_OSC> oscAmplitude = {};
			std::array<float, NUMBER_OF_OSC> oscHarmonic = {};
			float frequency = 0.0f;
			float amp = 0.0f;
		};

		std::array<SynthParams, 2> params = {};
		std::atomic<int> currentParamIndex{0};
		std::array<float, NUMBER_OF_OSC> _oscPhase = {};

		WaveTable waveTable;
		float _baseFrequency = 0.0f;
		float _z1 = 0.0f;


		void Process(float* buffPtr, int numberOfFrames) override {
			if (waveTable.table.size() == 0) return;
			int currentParams = currentParamIndex.load();

			for (int frame = 0; frame < numberOfFrames * ReadSettings().Channels; frame += ReadSettings().Channels) {
				float value = 0.0f;

				for (int osc = 0; osc < params[currentParams].numberOfActiveOsc; osc++) {
					value += params[currentParams].oscAmplitude[osc] * waveTable.Lookup(_oscPhase[osc]);

					_oscPhase[osc] += params[currentParams].oscInc[osc];
					if (_oscPhase[osc] >= 1.0f) _oscPhase[osc] -= 1.0f;
				}

				value *= params[currentParams].amp;

				for (int channel = 0; channel < ReadSettings().Channels; channel++) {
					buffPtr[frame + channel] += value + _z1;
				}

				_z1 = value;
			}
		}

		void Reset() override {
			waveTable.BuildTable();

			for (int p = 0; p <= 1; p++) {
				float osc_mag = 0;
				params[p].numberOfActiveOsc = 0;

				for (int i = 0; i < NUMBER_OF_OSC; i++) {
					params[p].oscHarmonic[i] = i+1;
					params[p].oscInc[i] = _baseFrequency * params[p].oscHarmonic[i];
					params[p].oscAmplitude[i] = 1.0f/(i+1);

					if (params[p].oscInc[i] < 0.5f) {
						params[p].numberOfActiveOsc++;
						osc_mag += params[p].oscAmplitude[i];
					}
				}

				params[p].amp = 1.0f / osc_mag;
				params[p].frequency = _baseFrequency;
			}
		}
	};

	class FMOscSoundChain : public SoundChainBase {
	public:
		FMOscSoundChain(float frequency = 440.0f) : _baseFrequency(frequency) {}
		~FMOscSoundChain() {}

		float GetFrequency() {return _baseFrequency;}
		void SetFrequency(float frequency) {
			// for (int i = 0; i < NUMBER_OF_OP; i++) {
			// 	_oscInc[i] = frequency * _oscRatio[i] / ReadSettings().SampleRate;
			// }
			// _baseFrequency = frequency;

			FMEvent event = {
				FMEventType::Frequency,
				0,
				0,
				frequency,
			};
			eventBuffer.Push(event);
		}

		void SetRatio(int op, float ratio) {
			if (op >= NUMBER_OF_OP) return;

			// _oscRatio[op] = ratio;
			// _oscInc[op] = _baseFrequency * ratio / ReadSettings().SampleRate;
			FMEvent event = {
				FMEventType::Ratio,
				op,
				0,
				ratio,
			};
			eventBuffer.Push(event);
		}

		void SetMod(int op, int mod, float level) {
			if (op >= NUMBER_OF_OP) return;

			// _oscFreqMod[op][mod] = level;
			FMEvent event = {
				FMEventType::Mod,
				op,
				mod,
				level,
			};
			eventBuffer.Push(event);
		}

		void SetMix(int op, float level) {
			if (op >= NUMBER_OF_OP) return;

			// _oscMix[op] = level;
			FMEvent event = {
				FMEventType::Mix,
				op,
				0,
				level,
			};
			eventBuffer.Push(event);
		}

	private:
		enum FMEventType {
			None,
			Frequency,
			Ratio,
			Mod,
			Mix
		};

		struct FMEvent {
			FMEventType eventType;
			int op;
			int mod;
			float value;
		};

		static const int NUMBER_OF_OP = 6;

		RingBuffer<FMEvent, 2> eventBuffer = {};

		std::array<float, NUMBER_OF_OP> _oscPhase = {};
		std::array<float, NUMBER_OF_OP> _oscInc = {};
		std::array<float, NUMBER_OF_OP> _oscMix = {};
		std::array<float, NUMBER_OF_OP> _oscRatio = {};
		std::array<std::array<float, NUMBER_OF_OP>, NUMBER_OF_OP> _oscFreqMod = {};
		std::array<float, NUMBER_OF_OP> _oscZ1 = {};

		WaveTable waveTable;
		float _baseFrequency = 0.0f;


		void Process(float* buffPtr, int numberOfFrames) override {
			if (waveTable.table.size() == 0) return;

			ProcessEvents();

			for (int frame = 0; frame < numberOfFrames * ReadSettings().Channels; frame += ReadSettings().Channels) {
				float value = 0.0f;
				std::array<float, NUMBER_OF_OP> phaseMod;

				for (int osc = 0; osc < NUMBER_OF_OP; osc++) {
					phaseMod[osc] = 0.0f;
					for (int mod = 0; mod < NUMBER_OF_OP; mod++) {
						phaseMod[osc] += _oscFreqMod[osc][mod] * _oscZ1[mod];
					}
				}

				for (int osc = 0; osc < NUMBER_OF_OP; osc++) {
					float sample = waveTable.Lookup(_oscPhase[osc] + phaseMod[osc]);
					value += _oscMix[osc] * sample;
					_oscZ1[osc] = sample;

					_oscPhase[osc] += _oscInc[osc];
					if (_oscPhase[osc] >= 1.0f) _oscPhase[osc] -= 1.0f;
				}

				for (int channel = 0; channel < ReadSettings().Channels; channel++) {
					buffPtr[frame + channel] += value;
				}
			}
		}

		void Reset() override {
			waveTable.BuildTable();

			for (int i = 0; i < NUMBER_OF_OP; i++) {
				_oscRatio[i] = i+1;
				_oscInc[i] = _baseFrequency * _oscRatio[i] / ReadSettings().SampleRate;
				_oscMix[i] = i == 0 ? 0.5f : 0.0f;

				if (i < NUMBER_OF_OP - 1) {
					_oscFreqMod[i][i+1] = 0.3f;
				}
			}
		}

		void ProcessEvents() {
			while (!eventBuffer.isEmpty()) {

				FMEvent event = eventBuffer.Pop();

				switch (event.eventType) {
				case FMEventType::Frequency:
					for (int i = 0; i < NUMBER_OF_OP; i++) {
						_oscInc[i] = event.value * _oscRatio[i] / ReadSettings().SampleRate;
					}
					_baseFrequency = event.value;
					break;
				case FMEventType::Ratio:
					_oscRatio[event.op] = event.value;
					_oscInc[event.op] = _baseFrequency * _oscRatio[event.op] / ReadSettings().SampleRate;
					break;
				case FMEventType::Mod:
					_oscFreqMod[event.op][event.mod] = event.value;
					break;
				case FMEventType::Mix:
					_oscMix[event.op] = event.value;
					break;
				case FMEventType::None:
					break;
				}
			}
		}
	};

}
