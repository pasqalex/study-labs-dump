FROM mcr.microsoft.com/dotnet/sdk:10.0 AS builder

WORKDIR /src
COPY . .
RUN dotnet publish -c Release -o /out


FROM mcr.microsoft.com/dotnet/runtime:10.0

WORKDIR /app
COPY --from=builder /out .

ENTRYPOINT ["dotnet", "LabFirstAgain.dll"]
