#pragma once

#include "SoundChainPlatform.hpp"

#include <cstdlib>
#include <windows.h>
#include <audioclient.h>
#include <mmdeviceapi.h>
#include <wrl/client.h>

#pragma comment(lib, "Ole32.lib")

using Microsoft::WRL::ComPtr;

class WasapiSCP : public SoundChainPlatform {
public:

private:
	ComPtr<IAudioClient>       _audioClient;
	ComPtr<IAudioRenderClient> _renderClient;

	HANDLE _audioEvent = nullptr;
	HANDLE _renderThread = nullptr;

	UINT32 _bufferFrames = 0;
	std::atomic<bool> _running{false};

	static DWORD WINAPI RenderThreadEntry(void* context) {
		return static_cast<WasapiSCP*>(context)->RenderThread();
	}

	DWORD WasapiSCP::RenderThread() {
		while (_running.load()) {
			DWORD result = WaitForSingleObject(_audioEvent, INFINITE);
			if (result != WAIT_OBJECT_0) break;

			UINT32 padding = 0;
			if (FAILED(_audioClient->GetCurrentPadding(&padding))) break;

			UINT32 framesAvailable = _bufferFrames - padding;
			if (framesAvailable == 0) continue;

			BYTE* data = nullptr;
			if (FAILED(_renderClient->GetBuffer(framesAvailable, &data))) break;

			FillBuffer(reinterpret_cast<float*>(data), static_cast<int>(framesAvailable));
			if (FAILED(_renderClient->ReleaseBuffer(framesAvailable, 0))) break;
		}

		return 0;
	}

	void Setup() override {
		HRESULT err;

		err = CoInitializeEx(nullptr, COINIT_MULTITHREADED);
		if (FAILED(err) && err != RPC_E_CHANGED_MODE) {
			printf("COM initialization failed");
			exit(EXIT_FAILURE);
		}

		ComPtr<IMMDeviceEnumerator> enumerator;
		ComPtr<IMMDevice> device;

		err = CoCreateInstance( __uuidof(MMDeviceEnumerator), nullptr, CLSCTX_ALL, IID_PPV_ARGS(&enumerator));
		if (FAILED(err)) {
			printf("Failed to create device enumerator");
			exit(EXIT_FAILURE);
		}

		err = enumerator->GetDefaultAudioEndpoint( eRender, eConsole, &device);
		if (FAILED(err)) {
			printf("Failed to get default output device");
			exit(EXIT_FAILURE);
		}

		err = device->Activate( __uuidof(IAudioClient), CLSCTX_ALL, nullptr, reinterpret_cast<void**>(_audioClient.GetAddressOf()));
		if (FAILED(err)) {
			printf("Failed to activate audio client");
			exit(EXIT_FAILURE);
		}

		WAVEFORMATEX format{};
		format.wFormatTag = WAVE_FORMAT_IEEE_FLOAT;
		format.nChannels = static_cast<WORD>(GetSoundChainSettings().Channels);
		format.nSamplesPerSec = static_cast<DWORD>(GetSoundChainSettings().SampleRate);
		format.wBitsPerSample = 32;
		format.nBlockAlign = format.nChannels * sizeof(float);
		format.nAvgBytesPerSec = format.nSamplesPerSec * format.nBlockAlign;

		err = _audioClient->Initialize(AUDCLNT_SHAREMODE_SHARED, AUDCLNT_STREAMFLAGS_EVENTCALLBACK, 0, 0, &format, nullptr);
		if (FAILED(err)) {
			printf("Failed to initialize WASAPI stream");
			exit(EXIT_FAILURE);
		}

		err = _audioClient->GetService(__uuidof(IAudioRenderClient), reinterpret_cast<void**>(_renderClient.GetAddressOf()));
		if (FAILED(err)) {
			printf("Failed to obtain render client");
			exit(EXIT_FAILURE);
		}

		err = _audioClient->GetBufferSize(&_bufferFrames);
		if (FAILED(err)) {
			printf("Failed to obtain WASAPI buffer size");
			exit(EXIT_FAILURE);
		}

		_audioEvent = CreateEvent( nullptr, FALSE, FALSE, nullptr);
		if (!_audioEvent) {
			printf("Failed to create audio event");
			exit(EXIT_FAILURE);
		}

		err = _audioClient->SetEventHandle(_audioEvent);
		if (FAILED(err)) {
			printf("Failed to set WASAPI event handle");
			exit(EXIT_FAILURE);
		}
	}

	void Start() override {
		_running.store(true);

		_renderThread = CreateThread(nullptr, 0, &WasapiSCP::RenderThreadEntry, this, 0, nullptr);
		if (!_renderThread) {
			printf("Failed to create audio thread");
			exit(EXIT_FAILURE);
		}

		if (FAILED(_audioClient->Start())) {
			printf("Failed to start WASAPI stream");
			exit(EXIT_FAILURE);
		}
	}

	void End() override {
		_running.store(false);

		if (_audioEvent) SetEvent(_audioEvent);

		if (_renderThread) {
			WaitForSingleObject(_renderThread, INFINITE);
			CloseHandle(_renderThread);
			_renderThread = nullptr;
		}

		if (_audioClient) _audioClient->Stop();

		if (_audioEvent) {
			CloseHandle(_audioEvent);
			_audioEvent = nullptr;
		}

		_renderClient.Reset();
		_audioClient.Reset();

		CoUninitialize();
	}
};